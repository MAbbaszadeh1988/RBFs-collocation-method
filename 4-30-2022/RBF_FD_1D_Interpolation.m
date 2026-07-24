clear
close all
clc
format shorte
%%
N = 100;
h = 1/N;
x = 0:h:1;
x = x';
M = length(x);
%%
ns = 23; % Size of Influence domain
F = zeros(M,ns);
C = zeros(M,ns);
for i=1:M
    x_center = x(i);
    rd = sqrt((x_center-x(:)).^2);
    [rd,ix] = sort(rd);
    F(i,:) = rd(1:ns);
    C(i,:) = ix(1:ns);
end
%%
phi    = @(r,c) exp(-(c*r).^2);
drphi  = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
d2rphi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));

wx1 = zeros(M);  % The first-order differential matrix D_{x}
wx2 = zeros(M);  % The second-order differential matrix D_{xx}
%% Optimal shape parameter S. Sarrah
c = 5;
minK = 1e16;
maxK = 1e18;
dc = 0.001;
%%
for i=1:M
    pn = C(i,:);
    rx = x(pn)-x(pn)';
    D_local = sqrt(rx.^2);

    K = 1;
    while (K<minK || K>maxK)
        A_local=phi(D_local,c);
        [~,S,~]=svd(A_local);
        K = S(1,1)/S(ns,ns);
        if K<minK  
            c = c - dc;
        elseif K>maxK 
            c = c + dc;
        end
    end

    B1_local = drphi(sqrt((x(i)-x(pn)).^2),x(i)-x(pn),c);
    B2_local = d2rphi(sqrt((x(i)-x(pn)).^2),x(i)-x(pn),c);
    wx1(i,pn) = pinv(A_local)*B1_local;
    wx2(i,pn) = pinv(A_local)*B2_local;
end
%%
f = @(x) sin(x);
fx = @(x) cos(x);
fxx = @(x) -sin(x);

norm(wx1*f(x)-fx(x),inf)
norm(wx2*f(x)-fxx(x),inf)