clear 
clc
format shorte
close all
tic
%%
a = -1.5; b = 1.5;
N = 10;
c = 0.5;
h = (b-a)/N;
n=6;
[x,y]=meshgrid(a:h:b);
[teta,r]=cart2pol(x,y);
xd=[];yd=[];p=0;q=0;
t=0:pi/40:2*pi-pi/40;
xb=(1/(n^2)).*(1+2.*n+n.^2-(n+1).*cos(n.*t)).*cos(t);
yb=(1/(n^2)).*(1+2.*n+n.^2-(n+1).*cos(n.*t)).*sin(t);
for i=1:length(teta)*length(r)
    if (r(i)*cos(teta(i)))^2+(r(i)*sin(teta(i)))^2<((1/(n^2))*(1+2*n+n^2-(n+1)*cos(n.*teta(i))))^2-0.001
        p=p+1;
        xd(p)=r(i)*cos(teta(i));
        yd(p)=r(i)*sin(teta(i));
    end
end
%%
T = 1;
tau = 1e-3;
M = T/tau;
t = 0:tau:T;
%%
u = @(x,y,t) exp(t).*sin(x+y);
f = @(x,y,t) 3.*exp(t).*sin(x+y);
%%
X = [xb';xd'];
Y = [yb';yd'];
NB = length(xb);
NI = length(xd);
N = NB + NI;
%%
phi = @(r) exp(-c^2.*r.^2);
dxphi = @(r,rx) -2*c^2*rx.*exp(-c^2.*r.^2);
dxxphi = @(r,rx) (-2*c^2+4*c^4.*rx.^2).*exp(-c^2.*r.^2);

dyphi = @(r,ry) -2*c^2*ry.*exp(-c^2.*r.^2);
dyyphi = @(r,ry) (-2*c^2+4*c^4.*ry.^2).*exp(-c^2.*r.^2);

r = sqrt( (X-X').^2 + (Y-Y').^2 );
rx = X-X';
ry = Y-Y';
%% 
A = phi(r);
Axx = dxxphi(r,rx);
Ayy = dyyphi(r,ry);

A_left  =  A - 0.5*tau.*(Axx+Ayy);
A_right =  A + 0.5*tau.*(Axx+Ayy);
%%
A_left(1:NB,:) = A(1:NB,:);

lambda = pinv(A)*u(X,Y,0);
for n=1:M
    t = (n-0.5)*tau;
    F = A_right*lambda + tau*f(X,Y,t);
    F(1:NB) = u(X(1:NB),Y(1:NB),n*tau);
    lambda = pinv(A_left)*F;
end
norm(A*lambda-u(X,Y,T),inf)
% hold on
% plot3(X,Y,A*lambda,'o')
Uap = A*lambda;

%%
a = -1.5; b = 1.5;
N = 200;
c = 0.5;
h = (b-a)/N;
n=6;
[x,y]=meshgrid(a:h:b);
[teta,r]=cart2pol(x,y);
xd=[];yd=[];p=0;q=0;
for i=1:length(teta)
    for j=1:length(r)
        if (r(i,j)*cos(teta(i,j)))^2+(r(i,j)*sin(teta(i,j)))^2<=((1/(n^2))*(1+2*n+n^2-(n+1)*cos(n.*teta(i,j))))^2
            xd(i,j)=r(i,j)*cos(teta(i,j));
            yd(i,j)=r(i,j)*sin(teta(i,j));
        else
            xd(i,j) = NaN;
            yd(i,j) = NaN;
        end
    end
end

Z = griddata(X,Y,Uap,xd,yd);
pcolor(xd,yd,Z), shading interp










