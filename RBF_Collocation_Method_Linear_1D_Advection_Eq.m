%% RBF collocation method for 1D Heat Eq.
clear
clc
format shorte
%% Collocation Points
n = 100;
L0 = -5;  L1 = 5;
h = (L1-L0)/n;
x = L0:h:L1;
boundary = find( x==L0 | x==L1 );
Interior = find( x~=L0 & x~=L1 );
N = length(x);
%% Time discrete scheme
dt = 1e-2;
T = 1;
M = floor(T/dt);
%% Basis Functions and Required Matrices
Phi = @(r,c) exp(-(c*r).^2);
dxPhi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
dxxPhi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));

rx =  x'-x ;
D = abs(rx);  % abs(x-x');
c = 10;

A = Phi(D,c);
Ax = dxPhi(D,rx,c);
Axx = dxxPhi(D,rx,c);
%% Continuous initial Condition
x = x';
u = @(x,t) (sech(x-t)).^2;
f = @(x,t) 0.*x;
%% Discontinuous initial Condition
% x = x';
% u = @(x,t) 1.*(0<=x-t & x-t<=0.5) + 0.*(x-t>0.5 & x-t<=1);
% f = @(x,t) 0.*x;
% % plot(x,u(x,T))
%%
A_left  = A + 0.5*dt*Ax;
A_right = A - 0.5*dt*Ax;

A_left(1,:) = A(1,:);

lambda = pinv(A)*u(x,0);

for n = 1:M
    
    F = A_right*lambda + dt*f(x,(n-0.5)*dt);
    
%     F(1) = 1;   % Discontinuous
    
    F(1) = u(x(1),n*dt); % Continuous
    
    lambda = pinv(A_left)*F;
    
    U = A*lambda;
    
    plot(x,U,'o',x,u(x,n*dt)), drawnow
%     pause(0.1)
end

% U = A*lambda;
% 
% plot(x,U,x,u(x,T))