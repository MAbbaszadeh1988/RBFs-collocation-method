%% RBF collocation method & Method of Lines for 1D Heat Eq.
clear
clc
format shorte
%% Collocation Points
n = 20;
L0 = 0;  L1 = 1;
h = (L1-L0)/n;
x = L0:h:L1;
boundary = find( x==L0 | x==L1 );
Interior = find( x~=L0 & x~=L1 );
N = length(x);
%% Time discrete scheme
dt = 1e-5;
T = 1;
M = floor(T/dt);
%% Basis Functions and Required Matrices
%     Phi = @(r,c) exp(-(c*r).^2);
%     dxPhi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
%     dxxPhi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));
%

Phi = @(r,c) sqrt( 1 + (c*r).^2 );
drPhi = @(r,rx,c) c^2*rx./sqrt( 1 + (c*r).^2 );
d2rPhi = @(r,rx,c) (c^2*( 1 + (c*r).^2 ) -c^4*rx.^2)./( ( 1 + (c*r).^2 ).^1.5 );

rx =  x'-x ;
D = abs(rx);  % abs(x-x');
c = 10;

A = Phi(D,c);
Ax = drPhi(D,rx,c);
Axx = d2rPhi(D,rx,c);
%%
u = @(x,t) exp(-t).*(sin(pi.*x)+cos(pi.*x));
f = @(x,t) exp(-t).*(pi^2*cos(pi.*x) + pi^2*sin(pi.*x)) - exp(-t).*(cos(pi.*x) + sin(pi.*x));
%%
RHS = @(t,lambda)  A\(Axx*lambda + f(x',t));

lambda = pinv(A)*u(x',0);

for n = 1:M
    
    t = n*dt;
    k1 = dt*RHS( t,lambda );
    k2 = dt*RHS( t+0.5*dt,lambda + 0.5*k1 );
    k3 = dt*RHS( t+0.5*dt,lambda + 0.5*k2 );
    k4 = dt*RHS( t+dt,lambda + k3 );
    
    lambda = lambda + (k1+2*k2+2*k3+k4)./6;
    
    U = A*lambda;
    
    U(boundary) = u(x(boundary),t);
    
    lambda = pinv(A)*U;
    
end

U = A*lambda;

norm(U-u(x',T),inf)
