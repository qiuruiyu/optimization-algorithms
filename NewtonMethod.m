clear;clc;cla;
syms x1 x2 t
% 2 dimensinal problem 
f(x1, x2) = (x1-x2)^2 + (x1-1)^2;
dfx1 = diff(f, x1);
dfx2 = diff(f, x2);
df = [dfx1; dfx2];
H = [diff(dfx1, x1) diff(dfx1, x2); diff(dfx2, x1), diff(dfx2, x2)];

x = linspace(-10, 10);
y = linspace(-10, 10);
[X1, X2] = meshgrid(x, y);
Z = (X1-X2).^2 + (X1-1).^2;
figure(1);
contour(X1, X2, Z);
axis equal;  
hold on;

% timing start
tic;

e = 1e-3;
Q = eye(2);
beta = 1;
alpha = 10;
iters = 0;
x0 = [0, 0];
x_record = x0;
scatter(x0(1), x0(2));
dfk = df(x0(1), x0(2));
while norm(dfk) > e
    iters = iters + 1;
    xk = x_record(end, :);
    if det(H(xk(1), xk(2))) ~= 0
        xk1 = xk - df(xk(1), xk(2))' * inv(H(xk(1), xk(2)));
    else
        G = H(xk(1), xk(2)) + beta * alpha * Q;
        beta = beta * alpha;
        xk1 = xk - df(xk(1), xk(2))' * inv(G);
    end
    x_record = [x_record; xk1];
    scatter(xk(1), xk(2));
    plot([xk(1), xk1(1)], [xk(2), xk1(2)], 'r', 'LineWidth', 2);
    dfk = df(xk1(1), xk1(2));
end

% timing end 
toc;

disp("Best solution found.");
disp("x = [" + double(xk1(1)) + ", " + double(xk1(2)) + "]");
disp("Iteration: " + iters);
disp("Running Time: " + toc);
