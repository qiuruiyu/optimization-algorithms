% CD Method 
% To solve f(x) = 1/2 x^TQx + b^Tx + c

clear;clc;cla; 
syms x1 x2 t 
% 2 dimensional problem 
x = [x1; x2];
A = randn(2);
% A = [1 0; 0 2];
Q = A'*A;
% Q = [2 0; 0 4];
b = randn(2, 1);
% b = [1; 2];
c = randn();
f(x1, x2) = 0.5*x'*Q*x + b'*x + c;
df = [diff(f, x1); diff(f, x2)];

x = linspace(-1, 1);
y = linspace(-1, 1);
[X1, X2] = meshgrid(x, y);
Z = 0.5 * (Q(1, 1)*X1.^2 + (Q(1, 2)+Q(2, 1))*X1.*X2 + (Q(2, 2)*X2.^2)) + b(1) * X1 + b(2) * X2 + c;
figure(1);
contour(X1, X2, Z);
axis equal;
hold on;

% timing start 
tic;

e = 1e-3;
iters = 0;
x0 = [0; 0];
x_record = x0;
scatter(x0(1), x0(2));
dfk = df(x0(1), x0(2));
dk = -dfk;
while norm(dfk) > e
    iters = iters + 1;
    xk = x_record(:, end);
    tk = - dk' * df(xk(1), xk(2)) / (dk' * Q * dk);
    xk1 = xk + tk * dk;
    scatter(xk1(1), xk1(2));
    plot([xk(1), xk1(1)], [xk(2), xk1(2)], 'r', 'LineWidth', 2);
    x_record = [x_record, xk1];
    dfk1 = df(xk1(1), xk1(2));
    if norm(dfk1) < e
        break;
    end
    beta = norm(dfk1)^2 / norm(dfk)^2;
    dk = -dfk1 + beta * dk;
    dfk = dfk1;
end

% timing end
toc;

disp("Best solution found.");
disp("x = [" + double(xk1(1)) + ", " + double(xk1(2)) + "]");
disp("Iteration: " + iters);
disp("Running Time: " + toc);



