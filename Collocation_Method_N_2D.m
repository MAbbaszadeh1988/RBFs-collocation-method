clear
clc
format shorte
%% 
n = 20;
a = 0; b = 1; d = 1;
h = (b-a)/n;
x = a:h:b;
y = a:h:d;
[X,Y] = meshgrid(x,y);
X = X(:);   Y = Y(:);
N = length(X);

boundary_n = find(Y==d);
boundary_s = find(Y==a);

boundary_e = find(X==b & Y~=a & Y~=d);
boundary_w = find(X==a & Y~=a & Y~=d);

Interior = find( X~=a &  X~=b & Y~=a & Y~=d);

% plot(X(boundary_e),Y(boundary_e),'ro',X(boundary_w),Y(boundary_w),'bo')
% axis equal
%%
Phi = @(r,c) exp(-(c*r).^2);
drPhi = @(r,rh,c) -2*c^2*rh.*exp(-(c*r).^2);
drrPhi = @(r,rh,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rh.^2));

rx = X-X';
ry = Y-Y';

D = sqrt( rx.^2 + ry.^2 );

c = 4;
Int = Phi(D,c);

D1x = drPhi(D,rx,c);
D1y = drPhi(D,ry,c);

D2x = drrPhi(D,rx,c);
D2y = drrPhi(D,ry,c);
%%
u = @(x,y) cos(pi.*x).*cos(pi.*y);

ux = @(x,y) -pi.*sin(pi.*x).*cos(pi.*y);
uy = @(x,y) -pi.*cos(pi.*x).*sin(pi.*y);



f = @(x,y) -2*pi^2.*cos(pi.*x).*cos(pi.*y);
%%
H = D2x + D2y;

H(boundary_n,:) = D1y(boundary_n,:);
H(boundary_s,:) = -D1y(boundary_s,:); 

H(boundary_e,:) = D1x(boundary_e,:);
H(boundary_w,:) = -D1x(boundary_w,:); 

F = f(X,Y);

F(boundary_n) = uy(X(boundary_n),Y(boundary_n));
F(boundary_s) = -uy(X(boundary_s),Y(boundary_s));

F(boundary_e) = ux(X(boundary_e),Y(boundary_e));
F(boundary_w) = -ux(X(boundary_w),Y(boundary_w));


Lambda = H\F;

U = Int*Lambda;

norm(U-u(X,Y),inf)

% hold on 
% plot3(X,Y,U,'o')

% U = reshape(U,length(x),length(y));
% X = reshape(X,length(x),length(y));
% Y = reshape(Y,length(x),length(y));
% 
% mesh(X,Y,U)
% 
% axis equal