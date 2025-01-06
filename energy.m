load test_wo_ans.mat

% 獲取 mask 的尺寸
[m, n] = size(mask);

% 獲取 kernel 的尺寸
[km, kn] = size(kernel);

% 將 kernel 填充到與 mask 相同的尺寸
kernel_padded = padarray(kernel, [m-km, n-kn], 'post');


% 計算 FFT
FFT_mask = fft2(mask);
FFT_kernel = fft2(kernel_padded);

% 執行卷積（頻域相乘）
FFT_convolution = FFT_mask .* FFT_kernel;

% 進行反 FFT 得到空間域結果
E = ifft2(FFT_convolution);

% 取實部，因為計算中可能引入微小虛數部分
E = real(E);

% 視覺化結果
imagesc(E);
title('Energy Deposition (卷積結果)');
colormap('jet');
colorbar;

% 在 mask 周圍進行零填充，擴展邊界
%mask_padded = padarray(mask, [km-1, kn-1]);

% 重新執行上述卷積步驟
