%% RBF collocation method & Method of Nonlinear for 1D Burger's Eq.
clear
clc
format shorte
%% Collocation Points
n = 40;
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
% Phi = @(r,c) exp(-(c*r).^2);
% drPhi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
% d2rPhi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));
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
eps = 1e-2;
RHS = @(t,lambda)  A\(eps.*Axx*lambda - (A*lambda).*(Ax*lambda));

initcond = @(x) sin(2*pi*x) + 0.5*sin(pi*x);
plot(x,initcond(x));

lambda = pinv(A)*initcond(x');

for n = 1:M
    
    t = n*dt;
    k1 = dt*RHS( t,lambda );
    k2 = dt*RHS( t+0.5*dt,lambda + 0.5*k1 );
    k3 = dt*RHS( t+0.5*dt,lambda + 0.5*k2 );
    k4 = dt*RHS( t+dt,lambda + k3 );
    
    lambda = lambda + (k1+2*k2+2*k3+k4)./6;
    
    U = A*lambda;
    
    U(boundary) = 0;
    
    lambda = pinv(A)*U;
    
    if mod(n,50)==0
        plot(x,U), drawnow
    end
end