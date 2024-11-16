% 加載數據
data = load('test_wo_ans.mat');
kernel = data.kernel;  % 電子能量分佈 (PSF)
mask = data.mask;      % 光刻圖形 (mask)

% 初始化劑量分佈
[m, n] = size(mask);   % 光刻圖形的大小
dose = ones(m, n);     % 初始劑量分佈

% 參數
tolerance = 1e-5;      % 收斂的容忍度
max_iter = 100;        % 最大迭代次數
expansion_factor = 2;  % FMM計算的擴展因子

% 計算目標能量分佈
target_energy = mask;

% 加載數據
data = load('test_wo_ans.mat');
kernel = data.kernel;  % 電子能量分佈 (PSF)
mask = data.mask;      % 光刻圖形 (mask)

% 初始化劑量分佈
[m, n] = size(mask);   % 光刻圖形的大小
dose = ones(m, n);     % 初始劑量分佈

% 參數
tolerance = 1e-5;      % 收斂的容忍度
max_iter = 100;        % 最大迭代次數
learning_rate = 0.1;   % 劑量修正的學習率

% 正確處理kernel的大小
[kernel_m, kernel_n] = size(kernel);

% 開始迭代劑量校正
for iter = 1:max_iter
    % 計算曝光能量分佈：對劑量和kernel進行卷積運算
    exposure_energy = conv2(dose, kernel, 'same');
    
    % 計算與目標光刻圖形的差異
    error = mask - exposure_energy;
    
    % 劑量修正：根據錯誤來更新劑量分佈
    dose = dose + learning_rate * error;  
    
    % 確保劑量值非負
    dose(dose < 0) = 0;  
    
    % 計算最大誤差，檢查收斂條件
    max_error = max(abs(error(:)));
    fprintf('第 %d 次迭代的最大誤差: %f\n', iter, max_error);
    
    if max_error < tolerance
        fprintf('收斂於第 %d 次迭代。\n', iter);
        break;
    end
end

% 結果可視化
figure;
subplot(1, 3, 1);
imagesc(mask);
title('目標光刻圖形 (Mask)');
colorbar;

subplot(1, 3, 2);
imagesc(exposure_energy);
title('曝光能量分佈');
colorbar;

subplot(1, 3, 3);
imagesc(dose);
title('最終劑量分佈');
colorbar;
