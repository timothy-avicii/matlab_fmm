% 加載數據
load test_wo_ans.mat

global max_particles
max_particles = 64;

global node_num
node_num = 0;

% 初始分割設置
[mask_r, mask_c] = size(mask);
mask_r = mask_r - 1;
mask_c = mask_c - 1;
lam_mask = gcd(mask_r, mask_c);

init_row = mask_r / lam_mask;
init_col = mask_c / lam_mask;

% 初始化節點
node = struct('range', [], 'part', 0, 'level', 1, 'children', []);
node = initial_split(init_row, init_col, node, lam_mask, mask);

% 遞歸分割
node = recursive_split(node, mask);

% 初始分割函數
function node = initial_split(row, col, node, lam, mask)
    % 用GCD分割mask為初始節點
    global node_num
    for i = 1:row
        for j = 1:col
            range = [lam * (i - 1) + 1, lam * i + 1, lam * (j - 1) + 1, lam * j + 1];
            sub_mask = mask(range(1):range(2)-1, range(3):range(4)-1);
            part_count = nnz(sub_mask);
            
            % 建立節點
            node(i, j).range = range;
            node(i, j).part = part_count;
            node(i, j).level = 1;
            node(i, j).children = [];
            
            % 更新節點總數
            node_num = node_num + 1;
        end
    end
end

% 遞歸分割函數
function node = recursive_split(node, mask)
    % 檢查是否需要進一步分割
    global max_particles
    global node_num
    
    for i = 1:numel(node)
        % 如果節點粒子數超過限制，則分割
        if node(i).part > max_particles
            range = node(i).range;
            mid_x = floor((range(1) + range(2)) / 2);
            mid_y = floor((range(3) + range(4)) / 2);
            
            % 分割為四個子區域
            children = struct('range', [], 'part', 0, 'level', 0 , 'children', []);
            children(1).range = [range(1), mid_x, range(3), mid_y];
            children(2).range = [mid_x+1, range(2), range(3), mid_y];
            children(3).range = [range(1), mid_x, mid_y+1, range(4)];
            children(4).range = [mid_x+1, range(2), mid_y+1, range(4)];
            [children.level] = deal(node(i).level + 1);
            for k = 1:4
                sub_mask = mask(children(k).range(1):children(k).range(2)-1, ...
                                children(k).range(3):children(k).range(4)-1);
                children(k).part = nnz(sub_mask);
                node_num = node_num + 1;
            end
            
            % 遞歸處理子節點
            for k = 1:4
                if children(k).part > max_particles
                    children(k) = recursive_split(children(k), mask);
                end
            end
            
            % 更新當前節點的子節點
            node(i).children = children;
        end
    end
end




