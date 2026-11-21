clear
close all

% python env 指定python編譯器
pyenv('Version', 'C:\Users\atim9\anaconda3\envs\pytorch-env\pythonw.exe');
pyenv(ExecutionMode = 'OutOfProcess');

% 引入模組
py.importlib.import_module('numpy');
py.importlib.import_module('torch');

% 傳入的路徑參數 python檔案路徑 
py.os.chdir('C:\Users\atim9\Desktop\程式\新程式\.vscode\PY\');
disp(['Python current directory: ', char(py.os.getcwd())]);

% 將python程式引入 不用加附檔名
in_func = py.importlib.import_module('test');
py.importlib.reload(in_func); % 再次引入，確保是最新的檔案

data1 = [1 2 3 ; 4 5 6 ];
data1 = py.torch.from_numpy(py.numpy.array(data1)).to(py.torch.float64);

data2 = double(in_func.trial(data1).numpy());

%   test.py
%   def trial()
%       print(sum(a))
