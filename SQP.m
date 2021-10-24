clear;clc;cla;
close all;
syms x1 x2 v;
% f(x1, x2) = x1^2 + x2^2;
f(x1, x2) = -20*exp(-0.2*sqrt(0.5*(x1^2+x2^2))) - exp(0.5*cos(2*pi*x1)+0.5*cos(2*pi*x2)) + 20 + exp(1);
g(x1, x2) = -x2 + x1;
L(x1, x2, v) = f - v * g;

x = linspace(-0, 5, 100);
y = linspace(-0, 5, 100);
[X1, X2] = meshgrid(x, y);
Z = -20*exp(-0.2.*sqrt(0.5.*(X1.^2+X2.^2))) - exp(0.5*cos(2*pi*X1)+0.5*cos(2*pi*X2)) + 20 + exp(1);
figure(1);
contour(X1, X2, Z);
% surf(X1, X2, Z);
hold on;
axis equal;

xk = [3.5, 3.5];
scatter(xk(1), xk(2), 'filled');
vk = 1;
Lx(x1, x2, v) = [diff(L, x1); diff(L, x2)];
Lv(x1, x2, v) = diff(L, v);
Lxx(x1, x2, v) = [diff(diff(L, x1), x1) diff(diff(L, x1), x2); diff(diff(L, x2), x1) diff(diff(L, x2), x2)];
df(x1, x2) = [diff(f, x1); diff(f, x2)];
% ddf(x1, x2, v) = [diff(diff(f, x1), x1) diff(diff(f, x1), x2); diff(diff(f, x2), x1) diff(diff(f, x2), x2)];
dg(x1, x2) = [diff(g, x1); diff(g, x2)];
tol = 1e-3;
iters = 0;
while sqrt(sum(Lx(xk(1), xk(2), vk).^2))>tol || abs(Lv(xk(1), xk(2), vk))>tol
    iters = iters + 1;  
    G = [Lxx(xk(1), xk(2), vk) -dg(xk(1), xk(2)); -dg(xk(1), xk(2))' 0];
    tmp = double(G \ [-Lx(xk(1), xk(2), vk); g(xk(1), xk(2))]);
    xk = [xk(1)+tmp(1); xk(2)+tmp(2)];
    scatter(xk(1), xk(2), 'filled');
    vk = vk + tmp(3);
end


