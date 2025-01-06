tic; % 開始計時

global n
n = 4;

global max_particles
max_particles = 4*n^2;

augmented_mask = augment(mask);
function y = augment(mask)
    mask(2^13,2^13)=0; %2^13是最靠近且大於mask邊長大小的2的指數，目的是讓mask在分割時，不出現小數座標
    y = mask;
end
levelnum = 1;

[mask_quadtree]= initialize_quadtree(augmented_mask, levelnum, ...
                                            normalized_boundary(augmented_mask), [1,1]);

%% Initialization
function normalized_boundary = normalized_boundary(matrix)
    %座標化矩陣的邊界 產生的陣列   

    normalized_boundary = [-1 1 -1 1];
end   %如果cell可以保持原矩陣元素的順序和位置，就可以不用頂點， %暫時不考慮
      
%%
function [quadtree] = initialize_quadtree(particles, levelnum, boundary, coordinate)
    % particles: 可能包含粒子的節點(矩陣)
    % boundary: 節點的邊界 [xmin xmax ymin ymax]
    % max_particles: 每個節點中最多能容納的粒子數
    %%% cell 尚未能成功實現
    

    % 如果當前節點的粒子數量超過 max_particles，則分割區域
    count = nnz(particles(:)) ;
    % count = 0;
    % for i = 1: size(particles,1)
    %     for j = 1: size(particles,2)
    %         if particles(i,j) == 1
    %             count = count+1;
    %         end
    %     end
    % end
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
        
        %cell = {};%mat2cell(particles,size(particles,1),size(particles,2));
        new_particles_in_cell = mat2cell(particles, [fix(size(particles,1)/2), ...
                                        the_other_half_of_particles_rows ], ...
                                       [fix(size(particles,2)/2), ...
                                        the_other_half_of_particles_cols] );

        
        % new_cell = new_particles_in_cell;
        % cell{1,1} = new_cell(1,1);
        % cell{1,2} = new_cell(1,2);
        % cell{2,1} = new_cell(2,1);
        % cell{2,2} = new_cell(2,2);
        %似乎不能單獨把cell陣列的一個cell切割                

    
        quadtree = struct('children', {}, 'divided_particles', [], 'is_leaf_node', [], ...
                          'level', levelnum,'boundary', boundary, 'coordinate',coordinate) ;
        
        % % 分割成四個區域
        % [region1, region2, region3, region4] = split_boundary(boundary);

        % 根據位置將粒子分配到各個區域  
        %可能可以省略
        quadtree(1).divided_particles = new_particles_in_cell{1,1};
        quadtree(2).divided_particles = new_particles_in_cell{1,2};
        quadtree(3).divided_particles = new_particles_in_cell{2,1};
        quadtree(4).divided_particles = new_particles_in_cell{2,2};
        
        quadtree(1).is_leaf_node = 0;
        quadtree(2).is_leaf_node = 0;
        quadtree(3).is_leaf_node = 0;
        quadtree(4).is_leaf_node = 0;

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
        
        quadtree(1).boundary = new_boundary1;
        quadtree(2).boundary = new_boundary2;
        quadtree(3).boundary = new_boundary3;
        quadtree(4).boundary = new_boundary4;

        %處理座標
        mid_x_prime = size(particles,2) / 2 + 1;
        mid_y_prime = size(particles,1) / 2 + 1;

        new_cor1 = [coordinate(1), coordinate(2)];
        new_cor2 = [mid_x_prime, coordinate(2)];
        new_cor3 = [coordinate(1), mid_y_prime];
        new_cor4 = [mid_x_prime, mid_y_prime];

        quadtree(1).coordinate = new_cor1;
        quadtree(2).coordinate = new_cor2;
        quadtree(3).coordinate = new_cor3;
        quadtree(4).coordinate = new_cor4;

        % 遞迴地為每個子區域建立四叉樹
        quadtree(1).children = initialize_quadtree(quadtree(1).divided_particles, ...
            levelnum + 1 , new_boundary1, new_cor1);
        if isempty(quadtree(1).children)
            quadtree(1).is_leaf_node = 1;
        end

        quadtree(2).children = initialize_quadtree(quadtree(2).divided_particles, ...
            levelnum + 1 , new_boundary2, new_cor2);
        if isempty(quadtree(2).children)
            quadtree(2).is_leaf_node = 1;
        end

        quadtree(3).children = initialize_quadtree(quadtree(3).divided_particles, ...
            levelnum + 1 , new_boundary3, new_cor3);
        if isempty(quadtree(3).children)
            quadtree(3).is_leaf_node = 1;
        end

        quadtree(4).children = initialize_quadtree(quadtree(4).divided_particles, ...
            levelnum + 1 , new_boundary4, new_cor4);
        if isempty(quadtree(4).children)
            quadtree(4).is_leaf_node = 1;
        end

    else
        % 如果粒子數量在可接受範圍內，則不分割
        % 需要給個空集合才能在後續使用isempty!!
        % quadtree = struct('children', {}, 'divided_particles', particles, ...
        %                   'level', levelnum,'boundary', boundary) ;
        % quadtree(1).level = levelnum;
        % quadtree(1).boundary = boundary;
        % quadtree(1).divided_particles = particles;
        % quadtree(1).children = {};
        % quadtree(1).is_leaf_node = [];
        quadtree = {};
    end
end  



%%做U_list
divided_mask_clone = BFS_for_U_list_1(mask_quadtree, mask);

function y = BFS_for_U_list_1(root, mask)
    queue = {root};
    % divided_mask_clone = num2cell(zeros(2^13)) ;
    divided_mask_clone = cell(max(size(mask,1), size(mask,2) ) );

    while ~isempty(queue)
        % 從佇列中取出第一個節點
        currentNode = queue{1};
        queue(1) = []; % 將取出的節點移除

        % % 顯示當前節點
        % disp(currentNode);
       
        % 檢查並加入子節點到佇列中
        if ~isempty(currentNode(1).children)
            queue{end + 1} = currentNode(1).children;
        %%%%%%% else
        %%%%%%%     currentNode(1).is_leaf_node = 1;
        %%%%%%%沒成功是因為，MATLAB 中的參數傳遞是值傳遞!!!!!!!!!!!
        end
        %currentnode底下可能沒分支了，所以要判斷他是否有大於2的分支(即有分支)
        if  ~isempty(currentNode(2).children)
            queue{end + 1} = currentNode(2).children;
        end

        if  ~isempty(currentNode(3).children)
            queue{end + 1} = currentNode(3).children;
        end

        if  ~isempty(currentNode(4).children)
            queue{end + 1} = currentNode(4).children;
        end
        % 因為是值傳遞，所以下面幾句不能直接修改到mask quadtree，要另外想辦法
        % if isempty(currentNode(1).children)
        %     currentNode(1).is_leaf_node = 1;
        % end
        % 
        % if isempty(currentNode(2).children)
        %     currentNode(2).is_leaf_node = 1;
        % end
        % 
        % if isempty(currentNode(3).children)
        %     currentNode(3).is_leaf_node = 1;
        % end
        % 
        % if isempty(currentNode(4).children)
        %     currentNode(4).is_leaf_node = 1;
        % end
        
        if(nnz(currentNode(1).divided_particles) ~= 0) %*如果最後結果有問題，要考慮是否不移除0節點
            % i是x座標，j是y座標

            %%因為記憶體會不足，所以從實心全賦值，改成空心

            %記憶體還是不太夠，進一步改成只賦值頂點座標
            %正方形左上與右下
            divided_mask_clone{currentNode(1).coordinate(1), ...
                currentNode(1).coordinate(2)}...
                = currentNode(1);

             if (currentNode(1).coordinate(1)-1+ ...
                size(currentNode(1).divided_particles,2)) <= size(mask,2)...
                && (currentNode(1).coordinate(2)+ ...
                size(currentNode(1).divided_particles,1) - 1) <= size(mask,1)

                divided_mask_clone{currentNode(1).coordinate(1)-1+ ...
                size(currentNode(1).divided_particles,2), ...
                ...
                currentNode(1).coordinate(2)+ ...
                size(currentNode(1).divided_particles,1) - 1}...
                ...
                = currentNode(1);
            end
            %正方形左下與右上
            if currentNode(1).coordinate(2)-1+ ...
                size(currentNode(1).divided_particles,1) <= size(mask,1)

                divided_mask_clone{currentNode(1).coordinate(1) , ...
                ...
                currentNode(1).coordinate(2)-1+ ...
                size(currentNode(1).divided_particles,1)}...
                = currentNode(1);
            end
            if currentNode(1).coordinate(1)+...
                size(currentNode(1).divided_particles,2) - 1 <= size(mask,2)

            divided_mask_clone{currentNode(1).coordinate(1)+...
                size(currentNode(1).divided_particles,2) - 1 , ...
                ...
                currentNode(1).coordinate(2)}...
                = currentNode(1);
            end
        end

        if(nnz(currentNode(2).divided_particles) ~= 0)
            %正方形左上與右下
            divided_mask_clone{currentNode(2).coordinate(1), ...
                currentNode(2).coordinate(2)}...
                = currentNode(2);

            if currentNode(2).coordinate(1)-1+ ...
                size(currentNode(2).divided_particles,2) <= size(mask,2)...
                && currentNode(2).coordinate(2)+ ...
                size(currentNode(2).divided_particles,1) - 1 <= size(mask,1)

                divided_mask_clone{currentNode(2).coordinate(1)-1+ ...
                size(currentNode(2).divided_particles,2), ...
                ...
                currentNode(2).coordinate(2)+ ...
                size(currentNode(2).divided_particles,1) - 1}...
                ...
                = currentNode(2);
            end
            %正方形左下與右上
            if currentNode(2).coordinate(2)-1+ ...
                size(currentNode(2).divided_particles,1) <= size(mask,1)

                divided_mask_clone{currentNode(2).coordinate(1) , ...
                ...
                currentNode(2).coordinate(2)-1+ ...
                size(currentNode(2).divided_particles,1)}...
                = currentNode(2);
            end
            if currentNode(2).coordinate(1)+...
                size(currentNode(2).divided_particles,2) - 1 <= size(mask,2)

            divided_mask_clone{currentNode(2).coordinate(1)+...
                size(currentNode(2).divided_particles,2) - 1 , ...
                ...
                currentNode(2).coordinate(2)}...
                = currentNode(2);
            end
        end

        if(nnz(currentNode(3).divided_particles) ~= 0)
            %正方形左上與右下
            divided_mask_clone{currentNode(3).coordinate(1), ...
                currentNode(3).coordinate(2)}...
                = currentNode(3);

            if currentNode(3).coordinate(1)-1+ ...
                size(currentNode(3).divided_particles,2) <= size(mask,2)...
                && currentNode(3).coordinate(2)+ ...
                size(currentNode(3).divided_particles,1) - 1 <= size(mask,1)

                divided_mask_clone{currentNode(3).coordinate(1)-1+ ...
                size(currentNode(3).divided_particles,2), ...
                ...
                currentNode(3).coordinate(2)+ ...
                size(currentNode(3).divided_particles,1) - 1}...
                ...
                = currentNode(3);
            end
            %正方形左下與右上
            if currentNode(3).coordinate(2)-1+ ...
                size(currentNode(3).divided_particles,1) <= size(mask,1)

                divided_mask_clone{currentNode(3).coordinate(1) , ...
                ...
                currentNode(3).coordinate(2)-1+ ...
                size(currentNode(3).divided_particles,1)}...
                = currentNode(3);
            end
            if currentNode(3).coordinate(1)+...
                size(currentNode(3).divided_particles,2) - 1 <= size(mask,2)

            divided_mask_clone{currentNode(3).coordinate(1)+...
                size(currentNode(3).divided_particles,2) - 1 , ...
                ...
                currentNode(3).coordinate(2)}...
                = currentNode(3);
            end
        end

        if(nnz(currentNode(4).divided_particles) ~= 0)
            %正方形左上與右下
            divided_mask_clone{currentNode(4).coordinate(1), ...
                currentNode(4).coordinate(2)}...
                = currentNode(4);
            if currentNode(4).coordinate(1)-1+ ...
                size(currentNode(4).divided_particles,2) <= size(mask,2)...
                && currentNode(4).coordinate(2)+ ...
                size(currentNode(4).divided_particles,1) - 1 <= size(mask,1)

                divided_mask_clone{currentNode(4).coordinate(1)-1+ ...
                size(currentNode(4).divided_particles,2), ...
                ...
                currentNode(4).coordinate(2)+ ...
                size(currentNode(4).divided_particles,1) - 1}...
                ...
                = currentNode(4);
            end
            %正方形左下與右上
            if currentNode(4).coordinate(2)-1+ ...
                size(currentNode(4).divided_particles,1) <= size(mask,1)

                divided_mask_clone{currentNode(4).coordinate(1) , ...
                ...
                currentNode(4).coordinate(2)-1+ ...
                size(currentNode(4).divided_particles,1)}...
                = currentNode(4);
            end
            if currentNode(4).coordinate(1)+...
                size(currentNode(4).divided_particles,2) - 1 <= size(mask,2)

            divided_mask_clone{currentNode(4).coordinate(1)+...
                size(currentNode(4).divided_particles,2) - 1 , ...
                ...
                currentNode(4).coordinate(2)}...
                = currentNode(4);
            end
            
        end
    end
    y = divided_mask_clone;
end


[U_list, V_list] = U_List_step2(divided_mask_clone);
function [U_list, V_list] = U_List_step2(divided_mask_clone)
    % % 預分配一個足夠大的結構體陣列
    % max_nodes = 8192^2 /64;
    U_list = struct('node', {}, 'U_list_of_node', {});
    V_list = struct('node', {}, 'V_list_of_node', {});
    counter1 = 1;
    counter2 = 1;

    % 使用單位步長遍歷所有節點
    for i = 1:size(divided_mask_clone, 2)
        for j = 1:size(divided_mask_clone, 1)

            if ~isempty(divided_mask_clone{i,j})
                current_node = divided_mask_clone{i,j};

                neighbors = find_neighbors(divided_mask_clone, i, j, current_node);
                if ~isempty(neighbors)
                    U_list(counter1).node = [i j];
                    U_list(counter1).U_list_of_node = neighbors;
                    counter1 = counter1 + 1;
                end

                neighbors_of_neighbors = find_neighbors_of_neighbors(divided_mask_clone, i, j, current_node);
                if ~isempty(neighbors_of_neighbors)
                    V_list(counter2).node = [i j];
                    V_list(counter2).V_list_of_node = neighbors_of_neighbors;
                    counter2 = counter2 + 1;
                end

            end
        end
    end
    
    % 移除未使用的預分配空間
    U_list = U_list(1:counter1-1);
    V_list = V_list(1:counter1-1);
end

function neighbors = find_neighbors(divided_mask_clone, i, j, current_node)
    neighbors = {};
    counter = 1;
    
    % 定義8個相鄰方向的偏移量
    offsets = [
        -1 -1; 0 -1; 1 -1;
        -1  0;       1  0;
        -1  1; 0  1; 1  1
    ];
    
    node_size = size(current_node.divided_particles);
    
    for k = 1:size(offsets, 1)
        ni = i + offsets(k,1) * node_size(2);
        nj = j + offsets(k,2) * node_size(1);
        
        if ni >= 1 && ni <= size(divided_mask_clone,2) && ...
           nj >= 1 && nj <= size(divided_mask_clone,1) && ...
           ~isempty(divided_mask_clone{ni,nj})
            
            neighbor = divided_mask_clone{ni,nj};
            if neighbor.level == current_node.level
                neighbors{counter} = [ni nj];
                counter = counter + 1;
            end
        end
    end
end

function neighbors_of_neighbors = find_neighbors_of_neighbors(divided_mask_clone, i, j, current_node)
    neighbors_of_neighbors = {};
    counter = 1;
    
    % 定義相鄰的相鄰方向的偏移量
    offsets = [
  -3 -3; -2 -3; -1 -3; 0 -3; 1 -3; 2 -3; 3 -3;
  -3 -2; -2 -2; -1 -2; 0 -2; 1 -2; 2 -2; 3 -2;
  -3 -1; -2 -1;                    2 -1; 3 -1; 
  -3  0; -2  0;                    2  0; 3  0;
  -3  1; -2  1;                    2  1; 3  1; 
  -3  2; -2  2; -1  2; 0  2; 1  2; 2  2; 3  2;
  -3  3; -2  3; -1  3; 0  3; 1  3; 2  3; 3  3;
    ];
    
    node_size = size(current_node.divided_particles);
    
    for k = 1:size(offsets, 1)
        ni = i + offsets(k,1) * node_size(2);
        nj = j + offsets(k,2) * node_size(1);
        
        if ni >= 1 && ni <= size(divided_mask_clone,2) && ...
           nj >= 1 && nj <= size(divided_mask_clone,1) && ...
           ~isempty(divided_mask_clone{ni,nj})
            
            neighbor = divided_mask_clone{ni,nj};
            if neighbor.level == current_node.level
                neighbors_of_neighbors{counter} = [ni nj];
                counter = counter + 1;
            end
        end
    end
end
%no empty entries checked



% % 
% % 
% % Upward pass
% % function multipole_expansion = upward_pass(node)
% %     如果是葉節點，計算多極展開
% %     if isempty(node{})
% %         multipole_expansion = compute_multipole(node.particles);
% %     else
% %         否則，遞迴地計算每個子節點的多極展開，並合併結果
% %         multipole_expansion = 0;
% %         for i = 1:4
% %             child_expansion = upward_pass(node.children{i});
% %             multipole_expansion = merge_multipole(multipole_expansion, child_expansion);
% %         end
% %     end
% % end
% % 
% % 
% % function multipole = compute_multipole(particles)
% %     根據粒子的位置和劑量值計算多極展開
% %     此處可自訂 Chebyshev 插值或其他方法
% %     multipole = sum(particles(:, 3));  % 假設第三列是劑量值
% % end
% % 
% % 
% % function result = merge_multipole(parent_multipole, child_multipole)
% %     合併父節點和子節點的多極展開
% %     result = parent_multipole + child_multipole;
% % end
% % 
% % 
% % 
% 
% 
% 
% % %% Downward pass
% % function local_expansion = downward_pass(node, parent_local_expansion)
% %     % 計算當前節點的局部展開，從父節點向下傳遞
% %     if isempty(node.children)
% %         local_expansion = parent_local_expansion + compute_local(node.particles);
% %     else
% %         % 遞迴地對每個子節點進行局部展開的傳遞
% %         for i = 1:4
% %             node.children{i}.local_expansion = downward_pass(node.children{i}, parent_local_expansion);
% %         end
% %     end
% % end
% % 
% % function local = compute_local(particles)
% %     % 計算粒子的局部展開
% %     local = sum(particles(:, 3));  % 假設第三列是劑量值
% % end
% % 
% % 
% 
% 
% 
% 
% % %% Gather
% % function energy_distribution = gather(node)
% %     if isempty(node.children)
% %         % 如果是葉節點，返回能量分佈
% %         energy_distribution = compute_energy(node.particles);
% %     else
% %         % 否則，遞迴地收集每個子節點的能量分佈
% %         energy_distribution = 0;
% %         for i = 1:4
% %             energy_distribution = energy_distribution + gather(node.children{i});
% %         end
% %     end
% % end
% % 
% % function energy = compute_energy(particles)
% %     % 計算每個粒子的最終能量分佈
% %     % 假設能量與多極展開和局部展開的結果相關
% %     energy = sum(particles(:, 3));  % 假設第三列是劑量值
% % end
% 
elapsedTime = toc; % 結束計時並返回所花時間（秒）
fprintf('程式執行時間：%.2f 秒\n', elapsedTime);

