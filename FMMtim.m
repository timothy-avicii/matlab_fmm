load test_wo_ans.mat
global max_particles
max_particles = 64;
levelnum = 1;
cell = {};
[mask_quadtree, cell ]= initialize_quadtree(mask,levelnum,cell,cal_boundary(mask));

%% Initialization
function boundary = cal_boundary(matrix)
    % 自動給出矩陣的邊界陣列 [1 xmax 1 ymax]  
    boundary = [1 size(matrix,1) 1 size(matrix,2)];
end   %如果cell可以保持原矩陣元素的順序和位置，就可以不用頂點， %暫時不考慮
      %改以雙迴圈，其位置在-1 0 -1的格子列為ULIST
    
function [quadtree, cell] = initialize_quadtree(particles,levelnum,cell, boundary)
    % particles: 可能包含粒子的節點(矩陣)
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

U_list = struct('node', '', 'list', [] ); 
V_list = struct('node', '', 'list', [] );
pre_order_traversal(mask_quadtree);

function pre_order_traversal = pre_order_traversal(root_node)
    U_counter1 = 0;
    V_counter1 = 0;
    if root_node(1).children ~= {}
        U_counter1  = U_counter1  + 1;

        %%U_list: 與一個小方塊四周的任一邊界相鄰者 且層級相同者
        %如何對於一個節點，讀取或列出所有與其相同層級的節點?? done?
        for i = 1:4
            for j = 1:4
                if (root_node(i).level == root_node(1).level)...
                    && ...
                  ( ...
                    (root_node(i).boundary(1) == root_node(1).boundary(1)|| ...
                    root_node(i).boundary(2) == root_node(1).boundary(1)||...
                    root_node(i).boundary(1) == root_node(1).boundary(2)|| ...
                    root_node(i).boundary(2) == root_node(1).boundary(2)...
                    )&&...
                    (root_node(i).boundary(3) == root_node(1).boundary(3)|| ...
                    root_node(i).boundary(4) == root_node(1).boundary(3)|| ...
                    root_node(i).boundary(3) == root_node(1).boundary(4)|| ...
                    root_node(i).boundary(4) == root_node(1).boundary(4)...
                    )...
                  )

                    U_list(U_counter1).list= [root_node(i).children ];
                end
            end
        end
        U_list(U_counter1).node= sprintf("%s%c%lf%s", inputname(1), '(', U_counter1 ...
                                        , ").children")
        % U_list(U_counter1).list= [root_node(2).children root_node(3).children
        %                 root_node(4).children ];

        pre_order_traversal(root_node(1).children);

    end

    if root_node(2).children ~= {}
        U_counter1  = U_counter1  + 1;

        for i = 1:4
            for j = 1:4
                if (root_node(i).level == root_node(1).level)...
                    && ...
                  ( ...
                    (root_node(i).boundary(1) == root_node(2).boundary(1)|| ...
                    root_node(i).boundary(2) == root_node(2).boundary(1)||...
                    root_node(i).boundary(1) == root_node(2).boundary(2)|| ...
                    root_node(i).boundary(2) == root_node(2).boundary(2)...
                    )&&...
                    (root_node(i).boundary(3) == root_node(2).boundary(3)|| ...
                    root_node(i).boundary(4) == root_node(2).boundary(3)|| ...
                    root_node(i).boundary(3) == root_node(2).boundary(4)|| ...
                    root_node(i).boundary(4) == root_node(2).boundary(4)...
                    )...
                  )
                    
                    U_list(U_counter1).list= [root_node(i).children ];
                end
            end
        end


        
        U_list(U_counter1).node= sprintf("%s%c%lf%s", inputname(1), '(', U_counter1 ...
                                        , ").children")
        % U_list(U_counter1).list= [root_node(1).children root_node(3).children
        %                 root_node(4).children ];
        % 
        pre_order_traversal(root_node(2).children);

    end

    if root_node(3).children ~= {}

        for i = 1:4
            for j = 1:4
                if (root_node(i).level == root_node(1).level)...
                    && ...
                  ( ...
                    (root_node(i).boundary(1) == root_node(3).boundary(1)|| ...
                    root_node(i).boundary(2) == root_node(3).boundary(1)||...
                    root_node(i).boundary(1) == root_node(3).boundary(2)|| ...
                    root_node(i).boundary(2) == root_node(3).boundary(2)...
                    )&&...
                    (root_node(i).boundary(3) == root_node(3).boundary(3)|| ...
                    root_node(i).boundary(4) == root_node(3).boundary(3)|| ...
                    root_node(i).boundary(3) == root_node(3).boundary(4)|| ...
                    root_node(i).boundary(4) == root_node(3).boundary(4)...
                    )...
                  )

                    U_list(U_counter1).list= [root_node(i).children ];
                end
            end
        end
        U_counter1  = U_counter1  + 1;
        U_list(U_counter1).node= sprintf("%s%c%lf%s", inputname(1), '(', U_counter1 ...
                                        , ").children")
        % U_list(U_counter1).list= [root_node(1).children root_node(2).children
        %                 root_node(4).children ];

        pre_order_traversal(root_node(3).children);

    end

    if root_node(4).children ~= {}

        for i = 1:4
            for j = 1:4
                if (root_node(i).level == root_node(1).level)...
                    && ...
                  ( ...
                    (root_node(i).boundary(1) == root_node(4).boundary(1)|| ...
                    root_node(i).boundary(2) == root_node(4).boundary(1)||...
                    root_node(i).boundary(1) == root_node(4).boundary(2)|| ...
                    root_node(i).boundary(2) == root_node(4).boundary(2)...
                    )&&...
                    (root_node(i).boundary(3) == root_node(4).boundary(3)|| ...
                    root_node(i).boundary(4) == root_node(4).boundary(3)|| ...
                    root_node(i).boundary(3) == root_node(4).boundary(4)|| ...
                    root_node(i).boundary(4) == root_node(4).boundary(4)...
                    )...
                  )
                    
                    U_list(U_counter1).list= [root_node(i).children ];
                end
            end
        end
        U_counter1  = U_counter1  + 1;
        U_list(U_counter1).node= sprintf("%s%c%lf%s", inputname(1), '(', U_counter1 ...
                                        , ").children")
        % U_list(U_counter1).list= [root_node(1).children root_node(2).children
        %                 root_node(3).children ];

        pre_order_traversal(root_node(4).children);

    end
end



% f