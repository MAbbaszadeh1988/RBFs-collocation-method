clear
clc
n = 20;
L0 = 0;  L1 = 1;
h = (L1-L0)/n;
x = L0:h:L1;
N = length(x);
D = abs(x-x');
%%
f = @(x) sin(x);
Phi = @(r,c) sqrt( 1 + (c*r).^2 );
% Phi = @(r,c) 1./sqrt( 1 + (c*r).^2 );
% Phi = @(r,c) exp(-(c*r).^2);
c = 0.5;
A = Phi(D,c);
%%
m = 2;
P = ones(N,m);
xp = x';
for j=2:m
    P(:,j) = xp.^(j-1);
end
%%
M = [A P;P' zeros(m)];
F = [f(xp);zeros(m,1)];
lambda = M\F;
%%
[~,tau] = cheb(10,L0,L1);%  Chebyshev points
y = (tau(end:-1:1))';  %sort(randn(1,N));
D_new =  abs(y'-x);
A_new = Phi(D_new,c);
yp = y;
P = ones(length(y),m);
for j=2:m
    P(:,j) = yp.^(j-1);
end
M_new = [A_new P];
S = M_new*lambda;
%%
norm(S-f(y'),inf)