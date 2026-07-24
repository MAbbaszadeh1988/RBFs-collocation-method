clear
clc
format shorte
%% 
n = 10;
a = 0; b = 1; d = 2;
h = (b-a)/n;
x = a:h:b;
y = a:h:d;
[X,Y] = meshgrid(x,y);
X = X(:);   Y = Y(:);
N = length(X);
tic
boundary = find( X==a | X==b | Y==a | Y==d );
Interior = find( X~=a &  X~=b & Y~=a & Y~=d);
toc
plot(X(boundary),Y(boundary),'o',X(Interior),Y(Interior),'*')
axis equal
%%
Phi = @(r,c) exp(-(c*r).^2);
drPhi = @(r,rh,c) -2*c^2*rh.*exp(-(c*r).^2);
drrPhi = @(r,rh,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rh.^2));

rx = X-X';
ry = Y-Y';

D = sqrt( rx.^2 + ry.^2 );

c = 2;
Int = Phi(D,c);

D1x = drPhi(D,rx,c);
D1y = drPhi(D,ry,c);

D2x = drrPhi(D,rx,c);
D2y = drrPhi(D,ry,c);
%%
u = @(x,y) sin(x).*sin(y);
f = @(x,y) -2.*sin(x).*sin(y);
%%
H = D2x + D2y;
H(boundary,:) = Int(boundary,:); 

F = f(X,Y);
F(boundary) = u(X(boundary),Y(boundary));

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