%% RBF collocation method & Method of Lines for 1D Heat Eq.
clear
clc
format shorte
%% Collocation Points
%% 
n = 10;
a = 0; b = 1; d = 1;
h = (b-a)/n;
x = a:h:b;
y = a:h:d;
[X,Y] = meshgrid(x,y);
X = X(:);   Y = Y(:);
N = length(X);
tic
boundary = find( X==a | X==b | Y==a | Y==d );
Interior = find( X~=a &  X~=b & Y~=a & Y~=d);
%% Time discrete scheme
dt = 1e-4;
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

rx = X-X';
ry = Y-Y';

D = sqrt( rx.^2 + ry.^2 );

c = 2;
A = Phi(D,c);

Axx = d2rPhi(D,rx,c);
Ayy = d2rPhi(D,ry,c);

%%
u = @(x,y,t) exp(-t).*sin(pi.*x).*sin(pi.*y);
f = @(x,y,t) exp(-t).*(2*pi^2-1).*sin(pi.*x).*sin(pi.*y);
%%
RHS = @(t,lambda)  A\((Axx+Ayy)*lambda + f(X,Y,t));

lambda = pinv(A)*u(X,Y,0);

for n = 1:M
    
    t = n*dt;
    k1 = dt*RHS( t,lambda );
    k2 = dt*RHS( t+0.5*dt,lambda + 0.5*k1 );
    k3 = dt*RHS( t+0.5*dt,lambda + 0.5*k2 );
    k4 = dt*RHS( t+dt,lambda + k3 );
    
    lambda = lambda + (k1+2*k2+2*k3+k4)./6;
    
    U = A*lambda;
    
    U(boundary) = u(X(boundary),Y(boundary),t);
    
    lambda = pinv(A)*U;
    
end

U = A*lambda;

norm(U-u(X,Y,T),inf)
