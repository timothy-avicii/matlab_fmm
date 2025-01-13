% 初始化 Chebyshev 類
N = 4; % 切比雪夫節點數量
D = 1; % 維度
cheb = chebyshev(N, D);

% 顯示節點
disp('Chebyshev Nodes:');
disp(cheb.nodes);
%% 

% 測試相似度矩陣
a = [0.5; -0.5]; % 第一組點
b = [0.8; -0.2]; % 第二組點
S = cheb.similarity(a, b);
disp('Similarity Matrix:');
disp(S);
%% 

% 測試向下插值
down_coeffs = cheb.downwards_coeffs();
disp('Downwards Coefficients:');
disp(down_coeffs);

% 測試向上插值
up_coeffs = cheb.upwards_coeffs();
disp('Upwards Coefficients:');
disp(up_coeffs);
