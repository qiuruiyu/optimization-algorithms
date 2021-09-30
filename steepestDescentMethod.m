clear;clc;cla;
syms x1 x2 t
% only for 2 dimensional condition
f(x1, x2) = (x1-x2)^2 + (x1-1)^2;
dfx1 = diff(f, x1);
dfx2 = diff(f, x2);

x = linspace(-1, 1);
y = linspace(-1, 1);
[X1, X2] = meshgrid(x, y);
Z = (X1-X2).^2 + (X1-1).^2;
figure(1);
contour(X1, X2, Z);
axis equal;  
hold on;

% timing start
tic;

e = 1e-3;
iters = 0;
x0 = [0, 0];
x_record = x0;
scatter(x0(1), x0(2), 'b');
df = [dfx1(x0(1), x0(2)), dfx2(x0(1), x0(2))];
while norm(df) > e
    iters = iters + 1;
    xk = x_record(end, :);
    % get best direction 
    dk = -[dfx1(xk(1), xk(2)), dfx2(xk(1), xk(2))];
    dk = dk / norm(dk);
    % get best step length
    g(x1, x2, t) = f(x1+dk(1)*t, x2+dk(2)*t);  
    dgt = diff(g, t);
    tk = solve(dgt(xk(1), xk(2), t));
    % update the x position 
%     xk1 = xk + tk * dk;
    % modify the step to avoid vertical condition 
    xk1 = xk + 0.9 * tk * dk;
    x_record = [x_record; xk1];
    scatter(xk1(1), xk1(2), 'b');
    plot([xk(1), xk1(1)], [xk(2), xk1(2)], 'r', 'LineWidth', 2);
    df = [dfx1(xk1(1), xk1(2)), dfx2(xk1(1), xk1(2))];
end

% timing end 
toc;

disp("Best solution found.");
disp("x = [" + double(xk1(1)) + ", " + double(xk1(2)) + "]");
disp("Iteration: " + iters);
disp("Running Time: " + toc);


