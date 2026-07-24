clear
clc
format shorte
%%
N = 30;
a = -2; b = 2; d = 2;
h = (b-a)/N;
%%
n = 6;
[x,y]=meshgrid(a:h:b);
[teta,r]=cart2pol(x,y);
xd=[];
yd=[];
p=0;
t=0:pi/40:2*pi-pi/40;
xb=(1/(n^2)).*(1+2.*n+n.^2-(n+1).*cos(n.*t)).*cos(t);
yb=(1/(n^2)).*(1+2.*n+n.^2-(n+1).*cos(n.*t)).*sin(t);

for i=1:length(teta)*length(r)
    if (r(i)*cos(teta(i)))^2+(r(i)*sin(teta(i)))^2<((1/(n^2))*(1+2*n+n^2-(n+1)*cos(n.*teta(i))))^2-0.01
        p=p+1;
        xd(p)=r(i)*cos(teta(i));
        yd(p)=r(i)*sin(teta(i));
    end
end

X = [xb';xd'];
Y = [yb';yd'];

NB = length(xb);
NI = length(xd);
Nt = NB + NI;
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
u = @(x,y) cos(pi.*x).*cos(pi.*y);
f = @(x,y) -2*pi^2.*cos(pi.*x).*cos(pi.*y);
%%
H = D2x + D2y;
H(1:NB,:) = Int(1:NB,:);

F = f(X,Y);
F(1:NB) = u(X(1:NB),Y(1:NB));

Lambda = H\F;

U = Int*Lambda;

norm(U-u(X,Y),inf)

% hold on
% plot3(X,Y,U,'o')

%%
h = 0.005
[x,y]=meshgrid(a:h:b);
[teta,r]=cart2pol(x,y);
x=[];
y=[];

for i=1:length(teta)
    for j=1:length(teta)
        if (r(i,j)*cos(teta(i,j)))^2+(r(i,j)*sin(teta(i,j)))^2<=((1/(n^2))*(1+2*n+n^2-(n+1)*cos(n.*teta(i,j))))^2
            x(i,j)=r(i,j)*cos(teta(i,j));
            y(i,j)=r(i,j)*sin(teta(i,j));
        else
            x(i,j) = NaN;
            y(i,j) = NaN;
        end
    end
end
%%
Z = griddata(X,Y,U,x(:),y(:));
Ze = griddata(X,Y,U-u(X,Y),x(:),y(:));

U_m = reshape(Z,length(x),length(y));
E_m = reshape(Ze,length(x),length(y));

figure, pcolor(x,y,U_m), shading interp
figure, pcolor(x,y,abs(E_m)), shading interp
