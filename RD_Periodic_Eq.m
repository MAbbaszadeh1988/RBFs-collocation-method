%% RBF collocation method for 2D Diffusion-Reaction Eq.
clear
clc
format shorte
%%  Spatial direction
n = 30;
a = 0; b = 1; d = 1;
h = (b-a)/n;
x = a:h:b;
y = a:h:d;
[X,Y] = meshgrid(x,y);
X = X(:);   Y = Y(:);
N = length(X);
boundary = find( X==a | X==b | Y==a | Y==d );
Interior = find( X~=a &  X~=b & Y~=a & Y~=d);

west=find(X==a & Y~=a & Y~=d);
east=find(X==b & Y~=a & Y~=d);
south=find(Y==a & X~=a & X~=b);
north=find(Y==d & X~=a & X~=b);
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

c = 10;
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


A_left(boundary,:)=A(north,:);
A_left(south,:)=A(south,:);

% A_left(east,:)=Ax(east,:)-Ax(west,:);
% A_left(west,:)=A(west,:)-A(east,:);

lambda = pinv(A)*u(X,Y,0);

for n = 1:M
    
    n*dt

    F = A_right*lambda + dt*f(X,Y,(n-0.5)*dt);
    
    F([north;south]) = u(X([north;south]),Y([north;south]),n*dt);

    
    lambda = pinv(A_left)*F;

    U = A*lambda;

    U(east) = U(west);

    lambda = pinv(A)*U;
 
        
end
clc
U = A*lambda;

norm(U-u(X,Y,T),inf)

% hold on
% plot3(X,Y,U,'o')

% U = reshape(U,length(x),length(y));
% X = reshape(X,length(x),length(y));
% Y = reshape(Y,length(x),length(y));
% %
% mesh(X,Y,U)
%
% axis equal