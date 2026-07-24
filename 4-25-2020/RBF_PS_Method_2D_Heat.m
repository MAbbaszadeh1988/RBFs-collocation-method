%% RBF-PS method for 2D Heat Eq.
clear
clc
format shorte
%% Collocation Points
n = 40;
c = 8;
L0 = 0;
L1 = 1;
h = (L1-L0)/n;
x = L0:h:L1;
[X,Y] = meshgrid(x);
X = X(:);   Y = Y(:);
N = length(X);
boundary = find( X==L0 | X==L1 | Y==L0 | Y==L1 );
Interior = find( X~=L0 &  X~=L1 & Y~=L0 & Y~=L1);
%% Time discrete scheme
T = 1;
dt = 1e-4;
M = floor(T/dt);
%% Basis Functions and Required Matrices
Phi = @(r,c) exp(-(c*r).^2);
drPhi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
d2rPhi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));

rx = X-X';
ry = Y-Y';

D = sqrt( rx.^2 + ry.^2 );
A = Phi(D,c);

Axx = d2rPhi(D,rx,c);
Ayy = d2rPhi(D,ry,c);
%%
u = @(x,y,t) exp(-t).*sin(pi.*x).*sin(pi.*y);
f = @(x,y,t) exp(-t).*(2*pi^2-1).*sin(pi.*x).*sin(pi.*y);
%%
I = eye(N);
Ainv = pinv(A);
DL  = I - 0.5*dt*(Axx+Ayy)*Ainv;
DR  = I + 0.5*dt*(Axx+Ayy)*Ainv;

DL(boundary,:) = I(boundary,:);

DLinv = pinv(DL);
U = u(X,Y,0);
for n = 1:M
    F = DR*U + dt*f(X,Y,(n-0.5)*dt);

    F(boundary) = u(X(boundary),Y(boundary),n*dt);

    U = DLinv*F;
end
norm(U-u(X,Y,T),inf)
