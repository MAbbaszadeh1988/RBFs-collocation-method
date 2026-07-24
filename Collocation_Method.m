%----------------------------------------------------------------
%  We want to solve | u''(x)+u(x)=f(x) s.t 0<x<1 by collocation method
%                   |
%                   | u(0)=0,  u(1)=0
%----------------------------------------------------------------
clear
clc
format shorte
V = [4 6 8 10 12 14 16];
E = zeros(length(V),1);
for k=1:length(V)
    n = V(k);
    h = 1/n;
    x = 0:h:1;
    N = length(x);
    %% Basis function
    phi = @(x,j) x.^(j) - x.^(j+1);
    dxxphi = @(x,j) (j)*(j-1).*x.^(j-2) - j*(j+1).*x.^(j-1);
    %% Matrices
    A = zeros(n-1);
    M = zeros(n-1);
    for i=1:n-1
        for j=1:n-1
            A(i,j) = dxxphi(x(i+1),j);
            M(i,j) = phi(x(i+1),j);
        end
    end
    %%
    u = @(x) sin(x)/sin(1) -x;
    f = @(x) -x;
    C = (A+M)\f(x(2:end-1)');
    U = M*C;
    E(k) = norm(U-u(x(2:end-1)'),inf);
end
%% Plot
semilogy(V,E,'-rs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.7,0.7,0.7])
