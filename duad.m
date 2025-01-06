% 假設layout是一個N x N的方形矩陣
layout = [ 0 0 0 0 0 0 0 0 ;
           0 0 1 0 0 0 0 0 ;
           0 0 1 0 1 1 0 0 ;
           0 0 1 1 1 0 0 0 ;
           0 0 0 0 0 0 0 0 ;
           0 1 1 1 1 0 0 0 ;
           0 0 0 0 1 1 1 0 ;
           0 0 0 0 0 0 0 0 
                ]; % 生成0和1的矩陣

%% 

a = size(layout);
N = a(1);
% 標準化坐標
% 假設layout是一個N x N的方形矩陣

% 標準化坐標
[x, y] = meshgrid(linspace(-1, 1, N), linspace(-1, 1, N));

% 顯示原始layout
figure;
imagesc(layout);
%title('原始layout');
axis equal;
colorbar;

% 建立四叉樹結構
maxDepth = 4; % 最大深度
quadtree = createQuadtree(layout, x, y, maxDepth);

% 顯示四叉樹結構
figure;
hold on;
%title('四叉樹結構');
axis([-1 1 -1 1]);
grid on;
plotQuadtree(quadtree);
hold off;

% 函數：建立四叉樹
function node = createQuadtree(layout, x, y, depth)
    if depth == 0 || sum(layout(:)) <= 4
        node = struct('layout', layout, 'x', x, 'y', y, 'children', []);
        return;
    end
    
    % 劃分區域
    midX = (min(x(:)) + max(x(:))) / 2;
    midY = (min(y(:)) + max(y(:))) / 2;
    
    idx1 = x <= midX & y <= midY;
    idx2 = x > midX & y <= midY;
    idx3 = x <= midX & y > midY;
    idx4 = x > midX & y > midY;
    
    % 建立子節點
    node = struct('layout', [], 'x', [], 'y', [], 'children', []);
    if sum(layout(idx1)) > 4
        node.children{1} = createQuadtree(layout(idx1), x(idx1), y(idx1), depth - 1);
    end
    if sum(layout(idx2)) > 4
        node.children{2} = createQuadtree(layout(idx2), x(idx2), y(idx2), depth - 1);
    end
    if sum(layout(idx3)) > 4
        node.children{3} = createQuadtree(layout(idx3), x(idx3), y(idx3), depth - 1);
    end
    if sum(layout(idx4)) > 4
        node.children{4} = createQuadtree(layout(idx4), x(idx4), y(idx4), depth - 1);
    end
end

% 函數：顯示四叉樹
function plotQuadtree(node)
    if isempty(node.children)
        scatter(node.x(:), node.y(:), 'filled');
        return;
    end
    
    for i = 1:4
        if ~isempty(node.children{i})
            plotQuadtree(node.children{i});
        end
    end
end
