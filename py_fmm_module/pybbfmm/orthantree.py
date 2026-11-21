import torch
from . import sets, ragged
from aljpy import arrdict
import numpy as np

def underoccupied(source_idxs, target_idxs, terminal, capacity):
    # Figure out which boxes are oversubscribed and need to be further subdivided
    source_unique, source_counts = torch.unique(source_idxs, return_counts=True)
    target_unique, target_counts = torch.unique(target_idxs, return_counts=True)

    source_okay = torch.ones_like(terminal)
    source_okay[source_unique] = (source_counts <= capacity)
    target_okay = torch.ones_like(terminal)
    target_okay[target_unique] = (target_counts <= capacity)

    return source_okay & target_okay


def orthantree(scaled, capacity=8):
    """Constructs a :ref:`tree <presolve>` for the given :func:`~pybbfmm.scale`'d problem.

    This is a bit of a mess of a function, but long story short it starts with all the sources allocated to the root
    and repeatedly subdivides overfull boxes, constructing the various tree tensors as it goes.

    :param scaled: :func:`~pybbfmm.scale`'d problem.
    :param capacity: the max number of sources or targets per box.
    :return: A :ref:`tree <presolve>`.
    """
    D = scaled.sources.shape[1]
    
    #源點和目標點一起計算
    points = torch.cat([scaled.sources, scaled.targets])
    indices = points.new_zeros((len(points),), dtype=torch.long)
        
    tree = arrdict.arrdict(
        parents=indices.new_full((1,), -1), # parent節點
        depths=indices.new_zeros((1,)), # 樹的深度
        centers=points.new_zeros((1, D)), # 節點的中心座標
        terminal=indices.new_ones((1,), dtype=torch.bool), # 是否為葉節點
        children=indices.new_full((1,) + (2,)*D, -1), # child節點
        descent=indices.new_zeros((1, D))) # 表示於parent節點的分裂方向
    
    bases = 2**torch.flip(torch.arange(D, device=indices.device), (0,)) #將child節點的分割方向做壓縮
    subscript_offsets = sets.cartesian_product(torch.tensor([0, 1], device=indices.device), D)
    center_offsets = sets.cartesian_product(torch.tensor([-1, +1], device=indices.device), D)
    
    depthcounts = [torch.as_tensor([1], device=indices.device)] # 記錄每一層的節點數目
    
    depth = 0
    while True:
        # 紀錄源點和目標點是否被計算
        used, used_inv = torch.unique(indices, return_inverse=True)
        source_idxs, target_idxs = indices[:len(scaled.sources)], indices[-len(scaled.targets):]
        tree.terminal[used] = underoccupied(source_idxs, target_idxs, tree.terminal, capacity)[used] 
        # indices紀錄目前每個points所在的樹節點位址
        # used紀錄已經使用過的樹節點位址
        # termina紀錄樹節點位址是否是葉節點，如果不是，則需要分割。
        # used_inv紀錄目前每個樹節點有哪些points，相當於indices的相反
        used_is_active = ~tree.terminal[used] # 判斷哪些節點是非終止樹節點
        point_is_active = used_is_active[used_inv] # 判斷哪些points屬於非終止樹節點
        
        if not point_is_active.any():
            break
        # 需要進行分裂
        depth += 1

        active = used[used_is_active] # 選取需要分裂的樹節點
        active_inv = (used_is_active.cumsum(0) - used_is_active.long())[used_inv[point_is_active]] # 選取需要分裂的points
        first_child = len(tree.parents) + 2**D*torch.arange(len(active), device=active.device) # 計算新增加的樹節點位址
        point_offset = ((points[point_is_active] >= tree.centers[active][active_inv])*bases).sum(-1) 
        # 用以計算新增加的樹節點 center
        child_box = first_child[active_inv] + point_offset  
        indices[point_is_active] = child_box
        
        trailing_ones = (slice(None),) + (None,)*D
        tree.children[active] = first_child[trailing_ones] + (subscript_offsets*bases).sum(-1)
        
        centers = tree.centers[active][trailing_ones] + center_offsets.float()/2**depth
        descent = center_offsets[None].expand_as(centers)
        
        n_children = len(active)*2**D
        children = arrdict.arrdict(
            parents=active.repeat_interleave(2**D),
            depths=tree.depths.new_full((n_children,), depth),
            centers=centers.reshape(-1, D),
            descent=descent.reshape(-1, D),
            terminal=tree.terminal.new_ones((n_children,)),
            children=tree.children.new_full((n_children,) + (2,)*D, -1))
        tree = arrdict.cat([tree, children])

        depthcounts.append(n_children)


    tree['id'] = torch.arange(len(tree.parents), device=points.device)

    indices = arrdict.arrdict(
        sources=indices[:len(scaled.sources)],
        targets=indices[-len(scaled.targets):])

    depths = ragged.Ragged(
        torch.arange(len(tree.id), device=points.device),
        torch.as_tensor(depthcounts, device=points.device))
    
    return tree, indices, depths


def child_boxes(tree, indices, descent):
    subscripts = (descent + 1)/2
    # unbind(-1) 依最後一維（各座標軸）拆開；對 (n, D) 與 (D,) 都和原本的 .T 結果相同，
    # 但不會觸發 PyTorch 對非二維 tensor 使用 .T 的棄用警告（未來版本會變成錯誤）
    return tree.children[(indices, *subscripts.long().unbind(-1))]


def neighbour_boxes(tree, indices, directions):
    """Finds the neighbour of `indices` in the tree in the given direction.

    There's a non-vectorized (and easier to understand) version of this function 
    `here <https://stackoverflow.com/questions/32412107/quadtree-find-neighbor/61211884#61211884>`_.
    """
    #TODO: This can be framed as a recursive scheme and then as a dynamic programming scheme. 
    # Should save a factor of log(n)
    indices = torch.as_tensor(indices, dtype=tree.parents.dtype, device=tree.parents.device)
    directions = torch.as_tensor(directions, dtype=tree.parents.dtype, device=tree.parents.device)
    directions = directions[None].repeat((len(indices), 1))
    assert len(directions) == len(indices), 'There should be as many directions as indices'
    
    current = indices.clone()
    alive = [torch.ones_like(indices, dtype=torch.bool)]
    neighbour_descents = [] # 父節點的相鄰節點的子節點偏移
    
        
    #raise SystemExit
    while alive[-1].any():
        live = alive[-1] & (directions != 0).any(-1) & (current >= 0)
        alive.append(live)
        
        descent = tree.descent[current]
        neighbour_descents.append(descent*(1 - 2*directions.abs()))
        
        directions = (descent + directions).div(2).long() 
        current[live] = tree.parents[current[live]]
        
    for descent, live in zip(neighbour_descents[::-1], alive[::-1]):
        
        internal = ~tree.terminal[current] & (current >= 0) & live
           
        current[internal] = child_boxes(tree, current[internal], descent[internal])
  
    return current
'''
# 這個是作者在stack_overflow提供的另一個版本的neighbour_boxes
def neighbour_boxes(tree, node, direction):
    """Finds the neighbour of a node in an adaptive quadtree or it's D-dimensional
    generalization (an orthantree?).

    Be very careful with indexing when implementing this. Because it indexes  
    physical space with arrays, sometimes you need to think of coords in terms
    of Euclidean coords, and other times in terms of array coords.
    
    Args:
        tree: an object holding a bunch of node attributes, indexed by node:
            * `parent`, an (N,)-array of integers. The `n`th element gives the 
               index of `n`'s parent. The parent of the root is -1.
            * `children`, an ((N,) + (2,)*D)-array of integers. The `n`th slice 
               gives the indices of `n`'s children. The bottom-left child of node 3
               in a quadtree would be (3, 0, 0). 
            * `descent`, an (N, D)-array with elements from {-1, +1}. The `n`th
               row gives which direction node `n` lies in compared to its parent.
               For example, the left-bottom quadrant in a quadtree would be `(-1, -1)`.
            * `terminal`, an (N,)-array of booleans. The `n`th element is True
               if node `n` is a leaf.
        node: an integer, the index of the node you want to find the neighbour of.
        direction: a (D,)-array with elements from {-1, +1}
    
    Returns:
        An integer giving the index of the neighbouring node, or -1 if it doesn't 
        exist.
    """
    direction = np.asarray(direction)
    
    # Ascend to the common ancestor
    neighbour_descents = []
    while True:
        if (direction == 0).all() or node < 0:
            break
        node_descent = tree.descent[node]
        neighbour_descent = node_descent*(1 - 2*abs(direction))
        neighbour_descents.append(neighbour_descent)
        
        #direction = ((node_descent + direction)/2).astype(int)
        direction = ((node_descent + direction)/2).to(torch.int64)
        

        node = tree.parents[node]   
    
    # Descend to the neighbour 
    for neighbour_descent in neighbour_descents[::-1]:
        if tree.terminal[node] or node < 0:
            break
        node = tree.children[(node, *((neighbour_descent + 1)//2).unbind(-1))]
  
    return node
'''




def u_scheme(tree, neighbours):
    """Calculates the u-:ref:`scheme <presolve>`.
    """
    unique_neighbours = torch.sort(neighbours, 1, descending=True).values
    unique_neighbours[:, 1:][unique_neighbours[:, 1:] == unique_neighbours[:, :-1]] = -1

    pairs = torch.stack([tree.id[:, None].expand_as(neighbours), unique_neighbours], -1)
    pairs = pairs[(pairs >= 0).all(-1) & tree.terminal[pairs].all(-1)]

    partner_is_larger = tree.depths[pairs[:, 0]] > tree.depths[pairs[:, 1]]
    smaller_partners = torch.flip(pairs[partner_is_larger], (1,))
    pairs = torch.cat([pairs, smaller_partners])

    return ragged.from_pairs(pairs, len(tree.id), len(tree.id))

def v_scheme(tree, depths, directions, neighbours):
    """Calculates the v-:ref:`scheme <presolve>`.
    """
    D = tree.children.ndim-1
    nonzero_directions = (directions != 0).any(-1)
    descents = sets.flat_cartesian_product(torch.tensor([-1, +1], device=tree.id.device), D)

    # The v list is many times bigger than the other lists, so we'll loop rather than 
    # vectorize to preserve memory.
    result = []
    for d in nonzero_directions.nonzero().squeeze(1):
        colleagues = neighbours[tree.parents, d]
        
        for friend_descent in descents:
            friends = child_boxes(tree, colleagues, friend_descent) 

            for own_descent in descents:
                offset = (-own_descent + 4*directions[d] + friend_descent)/2
                if (offset.abs() <= 1).all(-1):
                    continue

                for depth in torch.arange(depths.domain, device=tree.id.device):
                    s = depths.slice(depth)
                    mask = (tree.descent[s] == own_descent).all(-1) & ~tree.terminal[colleagues[s]] & (colleagues[s] >= 0)
                    result.append(arrdict.arrdict(
                        boxes=tree.id[s][mask],
                        friends=friends[s][mask],
                        offset=offset,
                        depth=depth)) 

    return result

def w_pairs(tree, directions, neighbours):
    """Metadata needed for calculating w- and x-:ref:`scheme <presolve>`.
    """
    D = tree.children.ndim-1
    bs = tree.terminal.nonzero().squeeze(1)

    origins, colleagues, dirs = [], [], []
    for d, dns in zip(directions, neighbours.T):
        dns = dns[bs]
        is_colleague = (tree.depths[bs] == tree.depths[dns])
        valid = is_colleague & ~tree.terminal[dns]
        origins.append(bs[valid])
        colleagues.append(dns[valid])
        dirs.append(d[None].repeat_interleave(valid.sum(), 0))
    origins, colleagues, dirs = torch.cat(origins), torch.cat(colleagues), torch.cat(dirs, 0)

    pairs = [origins.new_empty((0, 2))]
    parents = colleagues
    while parents.nelement():
        friends = tree.children[parents].reshape(-1, 2**D)
        distant = (tree.descent[friends] == dirs[:, None, :]).any(-1)
        
        ps = torch.stack([origins[:, None].expand_as(friends), friends], -1)
        pairs.append(ps[distant])
        
        mask = ~distant & ~tree.terminal[friends]
        origins, parents = ps[mask].T
        dirs = dirs[:, None].repeat_interleave(2**D, 1)[mask]
    pairs = torch.cat(pairs)

    return pairs

def interaction_scheme(tree, depths):
    """Returns the datastructures needed to calculate the :ref:`interactions <presolve>` between boxes.
    
    The datastructures are pretty heterogeneous because, well, performance. They're set
    up so the data needed can be got at fast without blowing up the memory budget.

    :param tree: a :ref:`tree <presolve>`.
    :param depths: the :ref:`depths <presolve>` to go with the tree.
    :return: a :ref:`scheme <presolve>`.
    """
    
    D = tree.children.ndim-1
    directions = sets.flat_cartesian_product(torch.tensor([-1, 0, +1], device=tree.id.device), D)
    # 生成八個方向的方向向量，代表相鄰節點 
    neighbours = torch.stack([neighbour_boxes(tree, tree.id, d) for d in directions], -1)

    w = w_pairs(tree, directions, neighbours)
    
    return arrdict.arrdict(
        w=ragged.from_pairs(w, len(tree.id), len(tree.id)),
        x=ragged.from_pairs(w.flip((1,)), len(tree.id), len(tree.id)),
        u=u_scheme(tree, neighbours),
        v=v_scheme(tree, depths, directions, neighbours))