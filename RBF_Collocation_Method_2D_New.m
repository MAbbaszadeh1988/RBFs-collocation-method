clear
clc
format shorte
%% 
n = 40;
L0 = 0;
L1 = 1;
h = (L1-L0)/n;
x = L0:h:L1;
[X,Y] = meshgrid(x);
X = X(:);   Y = Y(:);
N = length(X);
tic
boundary = find( X==0 | X==1 | Y==0 | Y==1 );
Interior = find( X~=0 &  X~=1 & Y~=0 & Y~=1);
toc
plot(X(boundary),Y(boundary),'o',X(Interior),Y(Interior),'*')

X = [X(boundary);X(Interior)];
Y = [Y(boundary);Y(Interior)];

NB = length(X(boundary));
NI = length(X(Interior));
%%
Phi = @(r,c) exp(-(c*r).^2);
dxPhi = @(r,rh,c) -2*c^2*rh.*exp(-(c*r).^2);
dxxPhi = @(r,rh,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rh.^2));

rx = X-X';
ry = Y-Y';

D = sqrt( rx.^2 + ry.^2 );

c = 2;
Int = Phi(D,c);

D1x = dxPhi(D,rx,c);
D1y = dxPhi(D,ry,c);

D2x = dxxPhi(D,rx,c);
D2y = dxxPhi(D,ry,c);
%%
u = @(x,y) sin(x).*sin(y);
f = @(x,y) -2.*sin(x).*sin(y);
%%
H = D2x + D2y;
H(1:NB,:) = Int(1:NB,:); 

F = f(X,Y);
F(1:NB) = u(X(1:NB),Y(1:NB));

Lambda = H\F;

U = Int*Lambda;

norm(U-u(X,Y),inf)