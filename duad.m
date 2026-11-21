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
    % 只要象限內有點就要建立子節點 (點數 <= 4 時遞迴會直接回傳葉節點)
    % 原本用 > 4 判斷，會把點數 1~4 的象限整個丟掉
    node = struct('layout', [], 'x', [], 'y', [], 'children', {cell(1, 4)});
    idxs = {idx1, idx2, idx3, idx4};
    for k = 1:4
        if sum(layout(idxs{k})) > 0
            node.children{k} = createQuadtree(layout(idxs{k}), x(idxs{k}), y(idxs{k}), depth - 1);
        end
    end
end

% 函數：顯示四叉樹
function plotQuadtree(node)
    if isempty(node.children)
        scatter(node.x(:), node.y(:), 'filled');
        return;
    end
    
    for i = 1:numel(node.children)
        if ~isempty(node.children{i})
            plotQuadtree(node.children{i});
        end
    end
end
