clear all
clc
format shorte
%%
N = 100;
h = 1/N;
x = 0:h:1;
x = x';
M = length(x);
%%
ns = 5;
c = zeros(M,1);
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
phi = @(r,s) sqrt(r.^2+s.^2);
phix = @(r,rh,s) rh./sqrt(r.^2+s.^2);
phixx = @(r,rh,s) 1./sqrt(r.^2+s.^2)-(rh.^2)./(r.^2+s.^2).^(1.5);
%%
wx2 = zeros(M);
Intlocal = zeros(M);
c = 10;
minK = 1e16;
maxK = 1e18;
dc = 0.001;
for i=1:M
    pn = C(i,:);
    rx = x(pn)-(x(pn))';
    D = sqrt(rx.^2);

    K = 1; 
    while (K<minK || K>maxK)
        A_local=phi(D,c);
        [U,S,V]=svd(A_local);
        K = S(1,1)/S(ns,ns);
        if K<minK,  c = c - dc;
        elseif K>maxK, c = c + dc;
        end
    end
    %%
    m = 4;
    P_local = ones(ns,m);
    dP_local = zeros(1,m);
    d2P_local = zeros(1,m);
    xp = x(pn);
    for j=2:m
        P_local(:,j) = xp.^(j-1);
        dP_local(1,j) = (j-1).*x(i).^(j-2);
        if j~=2
            d2P_local(1,j) = (j-1)*(j-2).*x(i).^(j-3);
        end
    end
    %%
    M_local  = [A_local P_local;P_local' zeros(m)];
    D_local = phix(sqrt((x(i)-x(pn)).^2),x(i)-x(pn),c);

    Mx_local = [D_local;dP_local(1,:)'];
    ww = pinv(M_local)*Mx_local;
    wx2(i,pn) = ww(1:ns);
end
%%
f = @(x) sin(x);
fx = @(x) cos(x);

norm(wx2*f(x)-fx(x),inf)