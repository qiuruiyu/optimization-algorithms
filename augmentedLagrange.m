clear;clc;cla;
close all;

syms x1 x2;
v = 1.5;
c = 10000;
% v: lagrange parameter 
% c: augmented lagrange parameter 
% min f(x1, x2)
f(x1, x2) = (x1-3)^2*(4-x2); 
% f(x1, x2) = -20*exp(-0.2*sqrt(0.5*(x1^2+x2^2))) - exp(0.5*cos(2*pi*x1)+0.5*cos(2*pi*x2)) + 20 + exp(1);
% constraint
g(x1, x2) = x1 + x2 + 3;
% Lagrange Function 
L(x1, x2) = f(x1, x2) - v * g(x1, x2); 
% Augmented Lagrange Function 
M(x1, x2) = L(x1, x2) + c / 2 * (g(x1, x2)^2);
% Hessian Matrix of ALM
H(x1, x2) = [diff(diff(M, x1), x1), diff(diff(M, x1), x2);
                    diff(diff(M, x2), x1), diff(diff(M, x2), x2)];

% visualize 
x = linspace(-0, 5, 100);
y = linspace(-0, 5, 100);
[X1, X2] = meshgrid(x, y);
Z = (X1-3).^2.*(4-X2); 
% Z = -20*exp(-0.2.*sqrt(0.5.*(X1.^2+X2.^2))) - exp(0.5*cos(2*pi*X1)+0.5*cos(2*pi*X2)) + 20 + exp(1);
% contour(X1, X2, Z);
surf(X1, X2, Z);
hold on;
axis equal;

tol = 1e-3;
% correction factor 
gamma = 1.5;
iters = 0 ;
xk = [-2.1234; 2];
scatter(xk(1), xk(2), 'filled');
text(xk(1), xk(2), string(iters), 'FontSize', 14);
x_record = xk;

while abs(g(xk(1), xk(2))) > tol
    iters = iters + 1;
    M(x1, x2) = L(x1, x2) + c / 2 * (g(x1, x2)^2);
    x_tmp = solve(diff(M, x1), diff(M, x2), x1, x2);
    xk = [real(double(x_tmp.x1)); real(double(x_tmp.x2))];
    len = length(x_tmp.x1);
    if len > 1 % more than one stable point
        for i = 1 : len % choose the local minimum
            xt1 = real(double(x_tmp.x1));
            xt2 = real(double(x_tmp.x2));
            xt = [xt1(i); xt2(i)]; % xt = [xt1; xt2];
            [~, flag] = chol(H(xt(1), xt(2))); % judge positive definite 
            if flag ==0 % is positive definite
                xk = xt;
                break;
            end
        end
        if i == len
            error("No local minimum point found in all stable points.");
        end
    end
    % calculate the change rate of g(x) 
    if g(xk(1), xk(2)) / g(x_record(1, end), x_record(2, end)) >= tol
        c = gamma * c;
    end
    v = v - c * g(xk(1), xk(2));
    scatter(xk(1), xk(2), 'filled');
    text(xk(1), xk(2), string(iters), 'FontSize', 14);
    line([x_record(1, end), xk(1)], [x_record(2, end), xk(2)]);
    x_record = [x_record, xk];
end




