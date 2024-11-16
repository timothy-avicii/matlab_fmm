load test_wo_ans.mat


global max_particles
max_particles = 64;

global levelnum
levelnum = 1;
cell = {};  

[mask_r, mask_c ] = size(mask) ;
mask_r = mask_r - 1 ;
mask_c = mask_c - 1 ;

lam_mask = gcd( mask_r , mask_c );

init_row =  mask_r/lam_mask ;
init_col = mask_c/lam_mask;
init_node( 1 , 1 ) = struct( 'range' , [ ] , 'part' , 0 ) ;
global node_num
node_num = 0 ;
init_node = initial_split(init_row,init_col, init_node,lam_mask,mask) ;

function node=initial_split(row,col,node,lam,mask)
    %初始分割
    
    global node_num
    for i= 1 : row
        for j = 1 : col
            node( i , j ).range = [lam*(i - 1) lam*i lam*(j - 1)  lam*j ] ;
            count = mask( lam*(i - 1)+1 : lam*i+1 , lam*(j - 1)+1 : lam*j+1 ) ;
            node( i , j ).part = nnz(count) ;
            node_num = node_num + 1 ;
        end
    end

    %用最大公因數做初始分割 每個節點一定都是正方形
end
levelnum = levelnum + 1 ;


%% 
function recursive_split(node)
    %遞歸分割
    global levelnum
    global node_num
    for i= 1 : row
        for j = 1 : col
            node( i , j ).range = [lam*(i - 1) lam*i lam*(j - 1)  lam*j ] ;
            count = mask( lam*(i - 1)+1 : lam*i+1 , lam*(j - 1)+1 : lam*j+1 ) ;
            node( i , j ).part = nnz(count) ;
            node(i , j ).level = levelnum ;
            node_num = node_num + 1 ;
        end
    end

end


function [q1, q2, q3, q4] = splitMatrix(data)
    [rows, cols] = size(data);
    midRow = ceil(rows / 2);  % 向上取整
    midCol = ceil(cols / 2);  % 向上取整
    %避免分割的範圍不是正方形
    q1 = data(1:midRow, 1:midCol);         % 左上象限
    q2 = data(1:midRow, midCol+1:end);     % 右上象限
    q3 = data(midRow+1:end, 1:midCol);     % 左下象限
    q4 = data(midRow+1:end, midCol+1:end); % 右下象限
end





%%
U_list = struct('node', '', 'list', []); 
V_list = struct('node', '', 'list', []); 

% 開始進行先序遍歷
pre_order_traversal(mask_quadtree);

function pre_order_traversal(root_node)
    
    global U_list;
    global V_list;
    
    % 初始化計數器
    U_counter1 = 0;
    V_counter1 = 0;
    
    % 遍歷四叉樹的四個子節點
    for k = 1:4
        if ~isempty(root_node(k).children)
            
            U_counter1 = U_counter1 + 1 ;
            % 更新U_list (最近鄰域節點)
            U_list(U_counter1).node = sprintf('%s%c%d%s', inputname(1), '(', U_counter1, ').children');
            U_list(U_counter1).list = find_neighbor_nodes(root_node, k);  % 優化：提取邏輯
            
            % 更新V_list (遠場節點)
            V_counter1 = V_counter1 + 1;
            V_list(V_counter1).node = sprintf('%s%c%d%s', inputname(1), '(', V_counter1, ').children');
            V_list(V_counter1).list = find_far_field_nodes(root_node, k);  % 優化：提取邏輯
            
            % 遞迴遍歷 recursive
            pre_order_traversal(root_node(k).children);
        end
    end
end

function neighbor_nodes = find_neighbor_nodes(root_node, current_index)
    % 查找與當前節點相鄰的同層級節點 (最近鄰域)
    neighbor_nodes = [];
    current_node = root_node(current_index);
    for i = 1:4
        if i ~= current_index  % 跳過當前節點本身
            if root_node(i).level == current_node.level && is_adjacent(root_node(i), current_node)
                neighbor_nodes = [neighbor_nodes, root_node(i).children];
            end
        end
    end
end

function far_field_nodes = find_far_field_nodes(root_node, current_index)
    % 查找與當前節點遠場相互作用的節點
    far_field_nodes = [];
    current_node = root_node(current_index);
    for i = 1:4
        if i ~= current_index  % 跳過當前節點本身
            if root_node(i).level == current_node.level && ~is_adjacent(root_node(i), current_node)
                far_field_nodes = [far_field_nodes, root_node(i).children];
            end
        end
    end
end

function adjacent = is_adjacent(node1, node2)
    % 判斷兩個節點是否相鄰
    adjacent = (node1.boundary(1) == node2.boundary(1) || node1.boundary(2) == node2.boundary(1) || ...
                node1.boundary(1) == node2.boundary(2) || node1.boundary(2) == node2.boundary(2)) && ...
               (node1.boundary(3) == node2.boundary(3) || node1.boundary(4) == node2.boundary(3) || ...
                node1.boundary(3) == node2.boundary(4) || node1.boundary(4) == node2.boundary(4));
end







