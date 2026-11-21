clear
close all

load test_wo_ans.mat
load BSL.mat
% python env 指定python編譯器
pyenv('Version', 'C:\Users\atim9\anaconda3\envs\pytorch-env\pythonw.exe');
pyenv(ExecutionMode = 'OutOfProcess');

% 引入模組
py.importlib.import_module('numpy');
py.importlib.import_module('torch');

% 傳入的路徑參數 python檔案路徑 
py.os.chdir('C:\Users\atim9\Desktop\程式\matlab_fmm\py_fmm_module\');
disp(['Python current directory: ', char(py.os.getcwd())]);

% 將python程式引入 不用加附檔名
pybbfmm_func = py.importlib.import_module('matlab_bench');
py.importlib.reload(pybbfmm_func); % 再次引入，確保是最新的檔案

% 將mask即為輸入資料
layout = py.torch.from_numpy(py.numpy.array(mask)).to(py.torch.int32); 
tic
result_ana_matrix = pybbfmm_func.matlab_fmm(layout).numpy();
toc

ana_matrix = double(result_ana_matrix);
dig_matrix = ana_matrix >= threshold;

MAE_matrix = abs(ana_matrix - energy_image);
MER_matrix = MAE_matrix ./ energy_image;
MAE_value = max(max(MAE_matrix));
MER_value = max(max(MER_matrix));
%%
imagesc(dig_matrix)
%title('Energy Distribution');
colorbar;
colormap(hot);

energy_image = ana_matrix;
mask_pattern = mask;
save('FMMdata.mat', 'energy_image', 'mask_pattern','threshold');
