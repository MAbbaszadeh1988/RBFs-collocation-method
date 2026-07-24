%% RBF-PS method for 1D Heat Eq.
clear
clc
format shorte
%% Collocation Points
n = 40;
c = 5;

L0 = 0;  L1 = 1;
h = (L1-L0)/n;
x = L0:h:L1;
boundary = find( x==L0 | x==L1 );
Interior = find( x~=L0 & x~=L1 );
N = length(x);
%% Time discrete scheme
T = 1;
dt = 1e-4;
M = floor(T/dt);
%% Basis Functions and Required Matrices
Phi = @(r,c) exp(-(c*r).^2);
drPhi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
d2rPhi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));

rx =  x'-x ;
D = abs(rx);  % abs(x-x');

A = Phi(D,c);
Ax = drPhi(D,rx,c);
Axx = d2rPhi(D,rx,c);
%%
u = @(x,t) exp(-t).*(sin(pi.*x)+cos(pi.*x));
f = @(x,t) exp(-t).*(pi^2*cos(pi.*x) + pi^2*sin(pi.*x)) - exp(-t).*(cos(pi.*x) + sin(pi.*x));

%%
I = eye(N);
Ainv = pinv(A);
DL  = I - 0.5*dt*Axx*Ainv;
DR  = I + 0.5*dt*Axx*Ainv;

DLinv = pinv(DL);
U = u(x',0);
for n = 1:M
    F = DR*U + dt*f(x',(n-0.5)*dt);
    U = DLinv*F;
    U(boundary) = u(x(boundary),n*dt);
end
norm(U-u(x',T),inf)
