%% RBF collocation method for 2D Diffusion-Reaction Eq.
clear
clc
format shorte
%%  Spatial direction
n = 20;
a = 0; b = 1; d = 1;
h = (b-a)/n;
x = a:h:b;
y = a:h:d;
[X,Y] = meshgrid(x,y);
X = X(:);   Y = Y(:);
N = length(X);
boundary = find( X==a | X==b | Y==a | Y==d );
Interior = find( X~=a &  X~=b & Y~=a & Y~=d);
%% Time discrete scheme
dt = 1e-2;
T = 1;
M = floor(T/dt);
%%
Phi = @(r,c) exp(-(c*r).^2);
drPhi = @(r,rh,c) -2*c^2*rh.*exp(-(c*r).^2);
drrPhi = @(r,rh,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rh.^2));

rx = X-X';
ry = Y-Y';

D = sqrt( rx.^2 + ry.^2 );

c = 2;
A = Phi(D,c);

Ax = drPhi(D,rx,c);
Ay = drPhi(D,ry,c);

Axx = drrPhi(D,rx,c);
Ayy = drrPhi(D,ry,c);
%%
u = @(x,y,t) exp(-t).*sin(pi.*x).*sin(pi.*y);
f = @(x,y,t) exp(-t).*sin(pi.*x).*sin(pi.*y).*( -2 +2*pi^2 );
%%
A_left  = (1-0.5*dt).*A - 0.5*dt*(Axx+Ayy);
A_right = (1+0.5*dt).*A + 0.5*dt*(Axx+Ayy);

A_left(boundary,:) = A(boundary,:);

lambda = pinv(A)*u(X,Y,0);

for n = 1:M
    
    n*dt

    F = A_right*lambda + dt*f(X,Y,(n-0.5)*dt);
    
    F(boundary) = u(X(boundary),Y(boundary),n*dt);
    
    lambda = pinv(A_left)*F;
    
%     Ua = A*lambda;
%     
%     Ua = reshape(Ua,length(x),length(y));
%     X1 = reshape(X,length(x),length(y));
%     Y1 = reshape(Y,length(x),length(y));
%     
%     mesh(X1,Y1,Ua), drawnow
    
    
end
clc
U = A*lambda;

norm(U-u(X,Y,T),inf)

% hold on
% plot3(X,Y,U,'o')

% U = reshape(U,length(x),length(y));
% X = reshape(X,length(x),length(y));
% Y = reshape(Y,length(x),length(y));
%
% mesh(X,Y,U)
%
% axis equal