clear;clc;cla;
close all;

syms x1 x2;
% 2 dimensional problem 
f(x1, x2) = 2*x1^2 -x2^2+ x2;
% constraint: g = x1 - x2 = 0;
g(x1, x2) = x1 - x2 - 1;
 
% visualize 
x = linspace(-3, 3, 100);
y = linspace(-3, 3, 100);
[X1, X2] = meshgrid(x, y);
Z = 2*X1.^2 - X2.^2 + X2;
contour(X1, X2, Z);
hold on;
axis equal;


gamma = 1;
sigma = 100;
% gamma = sigma * gamma
tol = 1e-6;
iters = 0;
xk = [-1; -1];
scatter(xk(1), xk(2), 'filled');
text(xk(1), xk(2), string(iters), 'FontSize', 14);
x_record = xk;
while abs(g(xk(1), xk(2))) > tol
    F(x1, x2) = f(x1, x2) + gamma * (g(x1, x2))^2;
    iters = iters + 1;
    x_tmp = solve(diff(F, x1)==0, diff(F, x2)==0, x1, x2);
    xk = [double(x_tmp.x1); double(x_tmp.x2)];
    scatter(xk(1), xk(2), 'filled');
    text(xk(1), xk(2), string(iters), 'FontSize', 14);
    line([x_record(1, end), xk(1)], [x_record(2, end), xk(2)]);
    x_record = [x_record, xk];
    gamma = sigma * gamma;
end

