%% RBF-PS method for 1D Heat Eq.
clear
clc
format shorte
%% Collocation Points
n = 50;
c = 0.5;
L0 = -5;  L1 = 5;
h = (L1-L0)/n;
x = L0:h:L1;
boundary = find( x==L0 | x==L1 );
Interior = find( x~=L0 & x~=L1 );
N = length(x);
%% Time discrete scheme
T = 0.5;
dt = 1e-4;
M = floor(T/dt);
%% Basis Functions and Required Matrices
Phi = @(r,c) exp(-(c*r).^2);
drPhi = @(r,rh,c) -2*c^2*rh.*exp(-(c*r).^2);
d2rPhi = @(r,rh,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rh.^2));
d3rPhi = @(r,rh,c) 4*c^4*rh.*exp(-(c*r).^2)* + 8*c^4*rh.*exp(-(c*r).^2) - 8*c^6*(rh.^3).*exp(-(c*r).^2);

rx =  x'-x ;
D = abs(rx);  % abs(x-x');

A = Phi(D,c);
Ax = drPhi(D,rx,c);
Axx = d2rPhi(D,rx,c);
Axxx = d3rPhi(D,rx,c);

%% Example 1
lambda = 0.0001;
alpha = 0.001;
beta = 0.001;
k = 0.4;
gama = 0;
mu1 = 2*k*sqrt(alpha+beta^2);
mu2 = 2*k^2*(alpha + beta^2 -beta*sqrt(alpha+beta^2));

u = @(x,t) lambda + mu1*tanh(k.*(x-lambda*t) + gama);
v = @(x,t) mu2 - mu2*(tanh(k.*(x-lambda*t) + gama).^2);

%%
I = eye(N);
Dx = Ax/A;
Dxx = Axx/A;
Dxxx = Axxx/A;

RHS_u = @(u,v) -Dx*v - u.*(Dx*u) + beta*Dxx*u;
RHS_v = @(u,v) -v.*(Dx*u) - u.*(Dx*v) + beta*Dxx*v - alpha*Dxxx*u;

U = u(x',0);
V = v(x',0);


for n = 1:M
    
    t = n*dt;

    k1 = RHS_u( U,V );
    l1 = RHS_v( U,V );

    k2 = RHS_u( U + 0.5*dt*k1 , V + 0.5*dt*l1 );
    l2 = RHS_v( U + 0.5*dt*k1 , V + 0.5*dt*l1 );

    k3 = RHS_u( U + 0.5*dt*k2 , V + 0.5*dt*l2 );
    l3 = RHS_v( U + 0.5*dt*k2 , V + 0.5*dt*l2 );

    k4 = RHS_u( U + dt*k3 , V + dt*l3 );
    l4 = RHS_v( U + dt*k3 , V + dt*l3 );
    
    U = U + dt*(k1+2.*k2+2.*k3+k4)./6;
    V = V + dt*(l1+2.*l2+2.*l3+l4)./6;
     
    
    U(boundary) = u(x(boundary),t);
    
    V(boundary) = v(x(boundary),t);
    
end
norm(U-u(x',T),inf)
figure, plot(x,U,'o',x,u(x',T))
figure, plot(x,V,'o',x,v(x',T))