clear
close all
n = 10 ;
 [a, b, c ] = generate2DChebyshevNodesUniformWeights(n);


function [X, Y, weights] = generate2DChebyshevNodesUniformWeights(n)
    % 生成二維切比雪夫節點和權重
    k = 0:n-1;
    nodes = cos((2*k + 1) * pi / (2*n));
    [X, Y] = meshgrid(nodes, nodes);
    weights = pi / n * ones(n, n); % 均勻權重
end
%% 
clear
close all

n = 10;

chebyshev_vs_cosine(n)

function chebyshev_vs_cosine(n)
    % 生成一維切比雪夫節點

    a = -1;       % 下限
    b = 1;        % 上限
    
    % 生成切比雪夫節點
    k = 0:n-1;
    p = (2*k + 1) * pi / (2*n)
    x_chebyshev = cos(p) % [-1, 1] 的節點
    x_mapped = ((b - a) / 2) * x_chebyshev + ((b + a) / 2); % 映射到 [a, b]
    
    % 生成cosine函數
    x_cosine = linspace(a, b, 500);  
    y_cosine = cos(pi * (x_cosine - a) / (b - a)); 
    
    
    figure;
    hold on;
    plot(x_cosine, y_cosine, 'b-', 'LineWidth', 2, 'DisplayName', 'Cosine Function');
    scatter(x_mapped, cos(pi * (x_mapped - a) / (b - a)), 100, 'r', 'filled', 'DisplayName', 'Chebyshev Nodes');
    hold off;
    
    
    title('Chebyshev Nodes vs. Cosine Function');
    xlabel('x');
    ylabel('f(x)');
    legend('show');
    grid on;
end
%% 
% 切比雪夫多项式插值和函数逼近
clc; clear; close all;


n = 20; % 切比雪夫多項式的階數
f = @(x) 1 ./ (1 + 25 * x.^2); % 待逼近函数（Runge 函数）
x_plot = linspace(-1, 1, 500); % 用于绘制原函数和插值多项式的点

% 生成切比雪夫節點
k = 1:n;
cheb_nodes = cos((2*k - 1) * pi / (2*n)); % 切比雪夫节点
f_nodes = f(cheb_nodes); % 函数值在切比雪夫节点上的取值

% 插值多项式系数（使用 Lagrange 插值或多项式拟合）
coeffs = polyfit(cheb_nodes, f_nodes, n-1); % 擬合多項式
f_interp = polyval(coeffs, x_plot); % 插值多项式的值


f_values = f(x_plot);


figure;

% 绘制原函数
plot(x_plot, f_values, 'b-', 'LineWidth', 2, 'DisplayName', 'Original Function f(x)');
hold on;

% 绘制插值多项式
plot(x_plot, f_interp, 'r--', 'LineWidth', 2, 'DisplayName', 'Chebyshev Interpolation');

% 绘制切比雪夫节点
scatter(cheb_nodes, f_nodes, 100, 'k', 'filled', 'DisplayName', 'Chebyshev Nodes');

% 图表设置
title('Chebyshev Polynomial Interpolation');
xlabel('x');
ylabel('f(x)');
legend('Location', 'Best');
grid on;
hold off;
%% 
clear 
close all
% 参数设置
n = 10; % 切比雪夫多项式的阶数
x = linspace(-1, 1, 500); % 用于绘图的点

% 切比雪夫多项式定义
Tn = cos(n * acos(x)); % 计算 T_n(x)

% 计算切比雪夫多项式的零点
k = 1:n; % 零点索引
zeros_chebyshev = cos((2*k - 1) * pi / (2*n)); % 零点公式

% 验证零点：在零点处，T_n(x) 应为 0
Tn_zeros = cos(n * acos(zeros_chebyshev)); % 计算零点处的 T_n(x)

% 绘图
figure;

% 绘制切比雪夫多项式曲线
plot(x, Tn, 'b-', 'LineWidth', 2, 'DisplayName', ['T', num2str(n), '(x)']);
hold on;

% 绘制零点
scatter(zeros_chebyshev, Tn_zeros, 100, 'r', 'filled', 'DisplayName', 'Zero Points');

% 设置图表属性
title(['Chebyshev Polynomial T', num2str(n), '(x) and Its Zero Points']);
xlabel('x');
ylabel(['T', num2str(n), '(x)']);
legend('Location', 'Best');
grid on;
hold on
%% 
% 二維矩陣的大小
M = 10; % 行數
N = 10; % 列數

% 生成原始二維矩陣
matrix = rand(M, N);

% 步驟1: 將矩陣座標化
[x, y] = meshgrid(1:N, 1:M); % 建立二維座標

% 步驟2: 轉換為 Chebyshev 節點
% Chebyshev 節點公式: x_i = cos((2i-1)*pi/(2n))，i = 1, 2, ..., n
% 對於 M 和 N 分別計算
chebyshev_x = cos((2*(1:N)-1)*pi/(2*N));
chebyshev_y = cos((2*(1:M)-1)*pi/(2*M));

% 步驟3: 映射矩陣資料到 Chebyshev 節點
% 生成新的網格 (Chebyshev 節點的網格)
[cheb_x, cheb_y] = meshgrid(chebyshev_x, chebyshev_y);

% 插值: 將原始矩陣映射到 Chebyshev 節點網格上
% 使用 interp2 進行插值
chebyshev_matrix = interp2(x, y, matrix, cheb_x, cheb_y, 'linear');

% 結果輸出
disp('原始矩陣:');
disp(matrix);

disp('Chebyshev 節點矩陣:');
disp(chebyshev_matrix);

% 可視化對比
figure;
subplot(1, 2, 1);
imagesc(matrix);
title('原始矩陣');
colorbar;

subplot(1, 2, 2);
imagesc(chebyshev_matrix);
title('Chebyshev 節點矩陣');
colorbar;
