load test_wo_ans.mat
global max_particles
max_particles = 64;
levelnum = 1;
cell = {};  
%cell 用來存葉節點

[mask_quadtree, cell ]= initialize_quadtree(mask,levelnum,cell,cal_boundary(mask));

%% Initialization
function boundary = cal_boundary(matrix)
    % 自動給出矩陣的邊界陣列 [1 xmax 1 ymax]  
    boundary = [1 size(matrix,1) 1 size(matrix,2)];
end   
%如果cell可以保持原矩陣元素的順序和位置，就可以不用頂點， %暫時不考慮
%改以雙迴圈，其位置在-1 0 -1的格子列為ULIST

%% 


function [quadtree, cell] = initialize_quadtree(particles,levelnum,cell, boundary)
    % particles: 可能包含粒子的節點(矩陣) mask切割後的
    % boundary: 節點的邊界 [xmin xmax ymin ymax]
    % max_particles: 每個節點中最多能容納的粒子數
 
    % 如果當前節點的粒子數量超過 max_particles，則分割區域
    count = nnz(particles(:)) ;

    for i = 1: size(particles,1)
        for j = 1: size(particles,2)
            if particles(i,j) == 1
                count = count+1;
            end
        end
    end
    
    global max_particles
    if count > max_particles %使用全域變數前也需要先宣告是global
        % 分割成四個區域
        
        if mod( size(particles,1), 2 ) == 1 %處理分隔點為小數的問題
            the_other_half_of_particles_rows = fix(size(particles,1)/2)+1;
        else
            the_other_half_of_particles_rows = size(particles,1)/2;
        end
        if mod( size(particles,2), 2 ) == 1
            the_other_half_of_particles_cols = fix(size(particles,2)/2)+1;
        else
            the_other_half_of_particles_cols = size(particles,2)/2;
        end
        
        cell = {};%mat2cell(particles,size(particles,1),size(particles,2));
        new_cell = mat2cell(particles, [fix(size(particles,1)/2), ...
                                        the_other_half_of_particles_rows ], ...
                                       [fix(size(particles,2)/2), ...
                                        the_other_half_of_particles_cols] );
        cell{1,1} = {new_cell{1,1}};
        cell{1,2} = {new_cell{1,2}};
        cell{2,1} = {new_cell{2,1}};
        cell{2,2} = {new_cell{2,2}};

        %cell是否能保有原本mask的全部元素?？                  
        %我想應該可以 因為本來就是由mask直接切割來的 
   
        quadtree = struct('children', {}, 'divided_particles', [], ...
                          'level', levelnum,'boundary', boundary) ;
        
        % % 分割成四個區域
        % [region1, region2, region3, region4] = split_boundary(boundary);

        % 根據位置將粒子分配到各個區域  
        %可能可以省略
        quadtree(1).divided_particles = new_cell{1,1};
        quadtree(2).divided_particles = new_cell{1,2};
        quadtree(3).divided_particles = new_cell{2,1};
        quadtree(4).divided_particles = new_cell{2,2};
        
        quadtree(1).level = levelnum;
        quadtree(2).level = levelnum;
        quadtree(3).level = levelnum;
        quadtree(4).level = levelnum;


        %處理邊界
        %平均分割
        mid_x = (boundary(1) + boundary(2)) / 2;
        mid_y = (boundary(3) + boundary(4)) / 2;

        new_boundary1 = [boundary(1), mid_x, boundary(3), mid_y];
        new_boundary2 = [mid_x, boundary(2), boundary(3), mid_y];
        new_boundary3 = [boundary(1), mid_x, mid_y, boundary(4)];
        new_boundary4 = [mid_x, boundary(2), mid_y, boundary(4)];
        
        quadtree(1).boundary = [boundary(1), mid_x, boundary(3), mid_y];
        quadtree(2).boundary = [mid_x, boundary(2), boundary(3), mid_y];
        quadtree(3).boundary = [boundary(1), mid_x, mid_y, boundary(4)];
        quadtree(4).boundary = [mid_x, boundary(2), mid_y, boundary(4)];


        % 遞迴地為每個子區域建立四叉樹
        quadtree(1).children = initialize_quadtree(quadtree(1).divided_particles, ...
            levelnum + 1,cell{1,1},new_boundary1);
        quadtree(2).children = initialize_quadtree(quadtree(2).divided_particles, ...
            levelnum + 1,cell{1,2},new_boundary2);
        quadtree(3).children = initialize_quadtree(quadtree(3).divided_particles, ...
            levelnum + 1,cell{2,1},new_boundary3);
        quadtree(4).children = initialize_quadtree(quadtree(4).divided_particles, ...
            levelnum + 1,cell{2,2},new_boundary4);
    else
        % 如果粒子數量在可接受範圍內，則不分割
        % 需要給個空集合才能在後續使用isempty!!
        quadtree = struct('children', {}, 'particles', particles);
    end

end  
% 問題:只有cell難以判斷葉節點(後續演算法需要用)，那麼可能還是要用struct

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



