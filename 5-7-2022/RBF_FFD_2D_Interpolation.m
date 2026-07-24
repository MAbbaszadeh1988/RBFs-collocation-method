clear
close all
clc
format shorte
%%
L0 = 0;   L1 = 1;
N = 80;
h = (L1-L0)/N;
x = L0:h:L1;
x = x';
[X,Y] = meshgrid(x);
X = X(:);   Y = Y(:);
Boundary = find(X==L0 | X==L1 | Y==L0 | Y==L1);
Interior = find(X~=L0 & X~=L1 & Y~=L0 & Y~=L1);
M = length(X);
%%
ns = 13; % Size of Influence domain
F = zeros(M,ns);
C = zeros(M,ns);
for i=1:M
    x_center = X(i);
    y_center = Y(i);
    rd = sqrt((x_center-X).^2+(y_center-Y).^2);
    [rd,ix] = sort(rd);
    F(i,:) = rd(1:ns);
    C(i,:) = ix(1:ns);
end
%%
phi    = @(r,c) exp(-(c*r).^2);
drphi  = @(r,rh,c) -2*c^2*rh.*exp(-(c*r).^2);
d2rphi = @(r,rh,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rh.^2));

wx1 = zeros(M);  % The first-order differential matrix D_{x}
wx2 = zeros(M);  % The second-order differential matrix D_{xx}
wy1 = zeros(M);  % The first-order differential matrix D_{y}
wy2 = zeros(M);  % The second-order differential matrix D_{yy}

%% Optimal shape parameter S. Sarrah
c = 0.5;
minK = 1e12;
maxK = 1e14;
dc = 0.1;
%%
for i=1:M
    pn = C(i,:);
    rx = X(pn)-X(pn)';
    ry = Y(pn)-Y(pn)';
    D_local = sqrt(rx.^2+ry.^2);

    K = 1;
    while (K<minK || K>maxK)
        A_local=phi(D_local,c);
        [~,S,~]=svd(A_local);
        K = S(1,1)/S(ns,ns)
        if K<minK  
            c = c - dc;
        elseif K>maxK 
            c = c + dc;
        end
    end
    %%
    I_local = sqrt( (X(i)-X(pn)).^2 + (Y(i)-Y(pn)).^2 );
    B1_local = drphi(I_local,X(i)-X(pn),c);
    B2_local = d2rphi(I_local,X(i)-X(pn),c);
    wx1(i,pn) = pinv(A_local)*B1_local;
    wx2(i,pn) = pinv(A_local)*B2_local;
    %%
    B1_local = drphi(I_local,Y(i)-Y(pn),c);
    B2_local = d2rphi(I_local,Y(i)-Y(pn),c);
    wy1(i,pn) = pinv(A_local)*B1_local;
    wy2(i,pn) = pinv(A_local)*B2_local;
end
%%
f = @(x,y) x.^3.*y.^3;
fx = @(x,y) (3.*x.^2) .*y.^3;
fxx = @(x,y) (6.*x).*(y.^3);

norm(wx1*f(X,Y)-fx(X,Y),inf)
norm(wx2*f(X,Y)-fxx(X,Y),inf)