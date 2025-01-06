% 定義切比雪夫節點
n = 4; % 切比雪夫多項式的階數
x = cos((2*(1:n)-1) * pi / (2*n));
y = cos((2*(1:n)-1) * pi / (2*n));

% 計算插值多項式 S_n(x, y)
[X, Y] = meshgrid(x, y);
Sn = cos(X) .* cos(Y);

% 顯示切比雪夫節點和插值多項式
figure;
scatter(X(:), Y(:), 'filled');
%title('切比雪夫節點');
xlabel('x');
ylabel('y');
grid on;

figure;
surf(X, Y, Sn);
%title('切比雪夫插值多項式 S_n(x, y)');
xlabel('x');
ylabel('y');
zlabel('S_n(x, y)');
grid on;
