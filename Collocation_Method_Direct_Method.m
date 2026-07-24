%----------------------------------------------------------------
%  We want to solve | u''(x)+u(x)=f(x) s.t 0<x<1 by collocation method
%                   |
%                   | u(0)=0,  u(1)=0
%----------------------------------------------------------------                 
clear
clc
format shorte
V = [4 6 8 10 12 16];
k = 1;
for n = V
h = 1/n;
x = 0:h:1;
N = length(x);
%% Basis function
phi = @(x,j) (1+x).^(j); %- x.^(j+1);
dxxphi = @(x,j) (j)*(j-1).*(1+x).^(j-2); %- j*(j+1).*x.^(j-1);
%% Matrices
A = zeros(N);
M = zeros(N);
I = eye(N);
for i=1:N
    for j=1:N
        A(i,j) = dxxphi(x(i),j);
        M(i,j) = phi(x(i),j);
    end
end
%%
u = @(x) sin(x)/sin(1) -x;
f = @(x) -x;
F = [u(0);f(x(2:end-1)');u(1)];
H = A+M;
H(1,:) = M(1,:);
H(end,:) = M(end,:);
C = H\F;
U = M*C;
E(k) = norm(U-u(x'),inf);
k = k+1;
end
%% Plot
semilogy(V,E,'-rs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.7,0.7,0.7])
