
clear 

% 初始化 Chebyshev 類
N = 3; % 切比雪夫節點數量
D = 3; % 維度

cheb = chebyshev(N, D);

% 顯示節點
disp('Chebyshev Nodes:');
disp(cheb.nodes);
% 
% % 相似度矩陣
% a = [0.5 0.1; -0.5 -0.1]; % 第一組點
% b = [0.8 0.5; -0.2 0.6]; % 第二組點
% S = cheb.similarity(a, b);
% disp('Similarity Matrix:');
% disp(S);

%向下插值
down_coeffs = cheb.downwards_coeffs();
disp('Downwards Coefficients:');
for i = 1 : numel(down_coeffs)
    down_coeffs(i)
end
%disp(down_coeffs);

% 向上插值
up_coeffs = cheb.upwards_coeffs();
disp('Upwards Coefficients:');
% for i = 1 : numel(up_coeffs)
%     up_coeffs(i)
% end
%disp(up_coeffs);
%%


clc;
clear;

% 初始化 Chebyshev 對象
N = 4; % 切比雪夫節點數量
D = 1; % 維度
cheb = chebyshev(N, D);

% 子節點數據
child_nodes = cheb.nodes; % Chebyshev 節點
v = [0.1 0.9  0.2 0.8]; % 子節點的初始數據

% Upward Pass: 計算父節點的權重
weights = cheb.anterpolate(child_nodes, v);
disp("Upward Pass (Weights):");
disp(weights)

% Downward Pass: 計算子節點的數據
interpolated_values = cheb.interpolate(child_nodes, weights);
disp("Downward Pass (Values):");
disp(interpolated_values)


% 誤差分析
error = norm(interpolated_values - v);
disp("Interpolation Error:");
disp(error);

%% 
clc;
clear;

% 初始化 Chebyshev 對象
N = 15;
D = 1;
cheb = chebyshev(N, D);

% 定義函數
g = @(x) 1 ./ (1 + 25 * x.^2);

% 節點與測試點
ns = cheb.nodes; % 切比雪夫節點
xs = linspace(-1, 1, 101)'; % 測試點

% 計算插值值
similarity_matrix = cheb.similarity(xs, ns); % 相似度矩陣
g_ns = g(ns); % 節點函數值
ghat = sum(similarity_matrix .* g_ns', 2); % 插值值

% 繪製比較
plot(xs, g(xs), 'b-', 'LineWidth', 1.5); % 真實函數
hold on;
plot(xs, ghat, 'r--', 'LineWidth', 1.5); % 插值近似
legend('Original Function', 'Interpolated Function');
title('Function Interpolation using Chebyshev Nodes');
xlabel('x');
ylabel('g(x)');
grid on;
%% 
% 初始化參數
N = 3; % Chebyshev nodes
D = 2; % Dimensions

% 建立 Chebyshev 物件
cheb = chebyshev(N, D);

% 定義測試函數
test_function = @(x) sum(x.^2, 2); % 測試函數為平方和

% 計算節點及其上的函數值
nodes = cheb.nodes; % Chebyshev 節點
function_values = test_function(nodes); % 函數值

% 隨機生成插值點
random_points = 2 * rand(5, D) - 1; % 在 [-1, 1]^D 隨機生成 5 點

% 使用 anterpolate 計算插值權重
weights = cheb.anterpolate(nodes, function_values);

% 使用 interpolate 計算插值點的值
interpolated_values = cheb.interpolate(random_points, weights);

% 計算真實值進行比較
true_values = test_function(random_points);

% 打印結果
disp('節點:');
disp(nodes);
disp('節點上的函數值:');
disp(function_values);
disp('隨機插值點:');
disp(random_points);
disp('插值計算值:');
disp(interpolated_values);
disp('真實值:');
disp(true_values);
