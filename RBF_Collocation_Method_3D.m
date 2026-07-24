clear
clc
format shorte
%% 
n = 10;
L0 = 0;
L1 = 1;
h = (L1-L0)/n;
x = L0:h:L1;
[X,Y,Z] = meshgrid(x);
X = X(:);   Y = Y(:);  Z = Z(:);
N = length(X);
tic
boundary = find( X==L0 | X==L1 | Y==0 | Y==1 | Z==0 | Z==1);
Interior = find( X~=L0 &  X~=L1 & Y~=0 & Y~=1 & Z~=0 & Z~=1);
%%
Phi = @(r,c) exp(-(c*r).^2);
drPhi = @(r,rh,c) -2*c^2*rh.*exp(-(c*r).^2);
drrPhi = @(r,rh,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rh.^2));

rx = X-X';
ry = Y-Y';
rz = Z-Z';

D = sqrt( rx.^2 + ry.^2 + rz.^2);

c = 2;
Int = Phi(D,c);

D1x = drPhi(D,rx,c);
D1y = drPhi(D,ry,c);
D1z = drPhi(D,rz,c);

D2x = drrPhi(D,rx,c);
D2y = drrPhi(D,ry,c);
D2z = drrPhi(D,rz,c);
%%
u = @(x,y,z) sin(x).*sin(y).*sin(z);
f = @(x,y,z) -3.*sin(x).*sin(y).*sin(z);
%%
H = D2x + D2y + D2z;
H(boundary,:) = Int(boundary,:); 

F = f(X,Y,Z);
F(boundary) = u(X(boundary),Y(boundary),Z(boundary));

Lambda = H\F;

U = Int*Lambda;

norm(U-u(X,Y,Z),inf)