%------------------------------------------------
% RBF-PS method for 1D Cubic Schrodinger  Eq.
% iu_t = -u_{xx} - 8abs(u)^2*u,  [-10,10]
% u(x,0) = sech(x)
%% ------------------------------------------------
clear
clc
format shorte
%% Collocation Points
n = 80;
c = 5;
L0 = -10;  L1 = 10;
h = (L1-L0)/n;
x = L0:h:L1;
boundary = find( x==L0 | x==L1 );
Interior = find( x~=L0 & x~=L1 );
N = length(x);
%% Time discrete scheme
T = 5;
dt = 1e-4;
M = floor(T/dt);
%% Basis Functions and Required Matrices
Phi = @(r,c) exp(-(c*r).^2);
d2rPhi = @(r,rh,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rh.^2));

rx =  x'-x ;
D = abs(rx);  % abs(x-x');

A = Phi(D,c);
Axx = d2rPhi(D,rx,c);
%%
I = eye(N);
Dxx = Axx/A;

RHS = @(u) (-Dxx*u - 8.*(abs(u).^2).*u)./1i;

U = sech(x');
H = zeros(M,N);
H(1,:) = U;
for n = 1:M

    t = n*dt;

    k1 = RHS( U );

    k2 = RHS( U + 0.5*dt*k1 );

    k3 = RHS( U + 0.5*dt*k2 );

    k4 = RHS( U + dt*k3 );

    U = U + dt*(k1+2.*k2+2.*k3+k4)./6;

    U(boundary) = 0;

    H(n+1,:) = U;

    if mod(n,100)==0
        plot(x,abs(U),'-o'), drawnow
    end
end

[x,t] = meshgrid(x,0:dt:T);
mesh(x,t,abs(H))
view(-10,45)