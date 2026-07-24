% RBF collocation method
clear
clc
format shorte
%% Collocation Points
n = 80;
L0 = 0;  L1 = 1;
h = (L1-L0)/n;
x = L0:h:L1;
N = length(x);
%%
% L0 = 0;  L1 = 1;
% [~,tau] = cheb(10,L0,L1);%  Chebyshev points
% x = (tau(end:-1:1))';  %sort(randn(1,N));
% N = length(x);
%% Basis Functions and Required Matrices
Phi = @(r,c) exp(-(c*r).^2);
dxPhi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
dxxPhi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));

rx =  x'-x ;
D = abs(rx);  % abs(x-x');
c = 2;

Int = Phi(D,c);
D1 = dxPhi(D,rx,c);
D2 = dxxPhi(D,rx,c);
%%
u = @(x) sin(pi.*x);
f = @(x) (1+pi^2).*sin(pi*x) + pi*cos(pi*x);
%% 
H = -D2 + D1 + Int; 
H(1,:) = D1(1,:) + Int(1,:); 
H(end,:) = D1(end,:) + Int(end,:); 

F = f(x');
F(1) = pi;
F(end) = -pi;

Lambda = H\F;

U = Int*Lambda;

norm(U-u(x'),inf)
