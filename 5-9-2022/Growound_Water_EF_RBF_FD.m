clear all
clc
format shorte
%%
tic
T = 1;
dt = 1e-4;
Nt = floor(T/dt);
%%
L0 = 0;   L1 = 1;
N = 20;
h = (L1-L0)/N;
x = L0:h:L1;
x = x';
D = 1e-2;
v = 1;
[X,Y] = meshgrid(x);
x = X(:);   y = Y(:);
Boundary = find(x==L0 | x==L1 | y==L0 | y==L1);
Interior = find(x~=L0 & x~=L1 & y~=L0 & y~=L1);
[M,~] = size(x);
[MB,~] = size(Boundary);
%% 
ns = 5;
c = zeros(M,1);
F = zeros(M,ns);
C = zeros(M,ns);
for i=1:M
    x_center = x(i);
    y_center = y(i);
    rd = sqrt((x_center-x(:)).^2+(y_center-y(:)).^2);
    [rd,ix] = sort(rd);
    F(i,:) = rd(1:ns);
    C(i,:) = ix(1:ns);
    c(i) = rd(ns)*sqrt(ns)*rd(2)/(0.02*(sum(F(i,1:ns))));
end
wx2 = zeros(M);
wy2 = zeros(M);
wx1 = zeros(M);
wy1 = zeros(M);
for i=1:M
    pn = C(i,:);
    rx = x(pn)-x(pn)';
    ry = y(pn)-y(pn)';
    r = sqrt(rx.^2+ry.^2);
    
    A_local = phi(r,c(i));
    
    Bx = D2(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2),x(i)-x(pn),c(i));
    wx2(i,pn) = A_local\Bx;
    By = D2(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2),y(i)-y(pn),c(i));
    wy2(i,pn) = A_local\By;
    
    Bx = D1(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2),x(i)-x(pn),c(i));
    wx1(i,pn) = A_local\Bx;
    By = D1(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2),y(i)-y(pn),c(i));
    wy1(i,pn) = A_local\By;
end
%%
Dx = wx1;  Dy = wy1;
Dxx = wx2;  Dyy = wy2;
%%

% S  = 0.3.*x.*(1-x).*y.*(1-y);
% A  = 0.2.*x.*(1-x).*y.*(1-y);
% Ms = 0.5.*x.*(1-x).*y.*(1-y);

%% Initial Conditions
U  = exp(-(x.^2+y.^2)./0.01);
V = U;
W = U;

pcolor(X,Y,reshape(U,size(X))), shading interp
%% 
RHS_U  = @(U,V,W) -v.*(Dx*U+Dy*U)+D.*(Dxx*U+Dyy*U)-0.6.*W.*(U./(1+U)).*(V./(2+V));
RHS_V  = @(U,V,W) -v.*(Dx*V+Dy*V)+D.*(Dxx*V+Dyy*V)-0.1.*W.*(U./(1+U)).*(V./(2+V));
RHS_W = @(U,V,W) -v.*(Dx*W+Dy*W)+D.*(Dxx*W+Dyy*W)-0.8.*W.*(U./(1+U)).*(V./(2+V))-2.*W;
m=1;
for jj = 1 : Nt
    jj*dt
    %%
    k1 = dt*RHS_U(U,V,W);
    k2 = dt*RHS_U(U+k1/2,V,W);
    k3 = dt*RHS_U(U+k2/2,V,W);
    k4 = dt*RHS_U(U+k3,V,W);
    U = U + (1/6)*( k1 + 2*(k2 + k3) + k4 );
    U(Boundary) = 0;
    %%
    k1 = dt*RHS_V(U,V,W);
    k2 = dt*RHS_V(U,V+k1/2,W);
    k3 = dt*RHS_V(U,V+k2/2,W);
    k4 = dt*RHS_V(U,V+k3,W);
    V = V + (1/6)*( k1 + 2*(k2 + k3) + k4 );
    V(Boundary) = 0;
    %%
    k1 = dt*RHS_W(U,V,W);
    k2 = dt*RHS_W(U,V,W+k1/2);
    k3 = dt*RHS_W(U,V,W+k2/2);
    k4 = dt*RHS_W(U,V,W+k3);
    W = W + (1/6)*( k1 + 2*(k2 + k3) + k4 );
    W(Boundary) = 0;
    
    [X1,Y1] = meshgrid(0:1/64:1);
    [Xq,Yq,Vq] = griddata(x,y,U,X1,Y1);
pcolor(Xq,Yq,Vq), shading interp
drawnow
end
[X1,Y1] = meshgrid(0:1/64:1);
[Xq,Yq,Vq] = griddata(x,y,U,X1,Y1);
pcolor(Xq,Yq,Vq), shading interp
colorbar on
% hold on
% plot3(x,y,S,'ro')
% U(:,:) = reshape(S,n,n);
% mesh(XX,YY,U); axis([0 1 0 1]); drawnow


% xx=[0,1];
% yy=0:h:1;
% [xb1,yb1]=meshgrid(xx,yy);
% xxx=h:h:(N-1)*h;
% yyy=[0,1];
% [xb2,yb2]=meshgrid(xxx,yyy);
% p = halton(20^2,2);
% xh = [xb1(:);xb2(:);p(:,1)];
% yh = [yb1(:);yb2(:);p(:,2)];
% % 
% [Xq,Yq,Vq] = griddata(x,y,S,xh,yh);
% tri = delaunay(Xq,Yq);
% trisurf(tri,Xq,Yq,Vq)
