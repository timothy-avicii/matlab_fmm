%load test_wo_ans.mat

%count = nnz(mask)

% 定義函數
f = @(x) 1 ./ (1 + 25 * x.^2);

% 節點數
n = 5;

% Chebyshev 節點
x_cheby = cos((2*(0:n-1)+1)/(2*n) * pi);

% 節點上的函數值
y_cheby = f(x_cheby);

% 插值點
x_interp = linspace(-1, 1, 1000);

% 使用 Chebyshev 插值進行插值
p = polyfit(x_cheby, y_cheby, n-1);
y_interp = polyval(p, x_interp);

% 繪圖
figure;
plot(x_interp, f(x_interp), 'b', 'DisplayName', '原函數');
hold on;
plot(x_interp, y_interp, 'r--', 'DisplayName', 'Chebyshev 插值');
plot(x_cheby, y_cheby, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Chebyshev 節點');
legend;
title('Chebyshev 插值示例');
xlabel('x');
ylabel('f(x)');
hold off;
%% 

%%測試資料

U_list = struct('node', '', 'list', []); 
V_list = struct('node', '', 'list', []); 

% 開始進行先序遍歷
pre_order_traversal(mask_quadtree);

function pre_order_traversal(root_node)
    % 初始化計數器
    U_counter1 = 0;
    V_counter1 = 0;
    
    % 遍歷四叉樹的四個子節點
    for k = 1:4
        if ~isempty(root_node(k).children)
            U_counter1 = U_counter1 + 1;
            
            % 更新U_list (最近鄰域節點)
            U_list(U_counter1).node = sprintf('%s%c%lf%s', inputname(1), '(', U_counter1, ').children');
            U_list(U_counter1).list = find_neighbor_nodes(root_node, k);  % 優化：提取邏輯
            
            % 更新V_list (遠場節點)
            V_counter1 = V_counter1 + 1;
            V_list(V_counter1).node = sprintf('%s%c%lf%s', inputname(1), '(', V_counter1, ').children');
            V_list(V_counter1).list = find_far_field_nodes(root_node, k);  % 優化：提取邏輯
            
            % 遞迴遍歷
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



%% 

% 定義結構體數組
U_list(1) = struct('node', 'U1', 'list', [1, 2, 3]);
U_list(2) = struct('node', 'U2', 'list', [4, 5, 6]);
U_list(3) = struct('node', 'U3', 'list', [7, 8, 9]);

% 顯示結構體數組
disp(U_list);
