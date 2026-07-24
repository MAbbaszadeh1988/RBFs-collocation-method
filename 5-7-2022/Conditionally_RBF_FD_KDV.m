clear all;
close all;
clc;
dt = 1e-1;
Tfin = 50;
nmax = floor(Tfin/dt);
L1 = 100;
tspan = 0 : dt: Tfin;
n = 400;
x = linspace ( -L1, L1, n )';
ns = 101;
M = length(x);
%%
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
drphi = @(r,rh,s) rh./sqrt(r.^2+s.^2);
d3rphi = @(r,rh,s) (3*rh.^3)./(s^2 + r.^2).^(5/2) - (3.*rh)./(s^2 + r.^2).^(3/2);

W1x = zeros(M);  % The first-order differential matrix D_{x}
W3x = zeros(M);  % The second-order differential matrix D_{xx}
%% Optimal shape parameter S. Sarrah
c = 2;
minK = 1e8;
maxK = 1e10;
dc = 0.1;
%%
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
    m = 2;
    P = ones(ns,m);
    dP = zeros(1,m);
    d2P = zeros(1,m);
    d3P = zeros(1,m);
    xp = x(pn);
    for j=2:m
        P(:,j) = xp.^(j-1);
        dP(1,j) = (j-1).*x(i).^(j-2);
        if j~=2
            d2P(1,j) = (j-1)*(j-2).*x(i).^(j-3);
        end
        if j~=3
            d3P(1,j) = (j-1)*(j-2)*(j-3).*x(i).^(j-4);
        end
    end
    M_local  = [A_local P;P' zeros(m)];
    D_local = drphi(sqrt((x(i)-x(pn)).^2),x(i)-x(pn),c);
    D3_local = d3rphi(sqrt((x(i)-x(pn)).^2),x(i)-x(pn),c);

    Mx_local  = [D_local;dP(1,:)'];
    M3x_local = [D3_local;d3P(1,:)'];

    ww = pinv(M_local)*Mx_local;
    W1x(i,pn) = ww(1:ns);

    ww = pinv(M_local)*M3x_local;
    W3x(i,pn) = ww(1:ns);

    
end
%%
%Zero flux boundary conditions
W1x(1,:) = zeros(size(W1x(1,:))); W1x(end,:) = zeros(size(W1x(1,:)));
W3x(1,:) = zeros(size(W1x(1,:))); W3x(end,:) = zeros(size(W1x(1,:)));
A = sqrt(1)/sqrt(6); L = 1; x0 = 0;
U = @(x,t) 3*A^2*sech(A*L*(x - x0/L)/2 - A^3*t/2).^2; %Exact solution
RHS_u = @(t,u) -u.*(W1x*u) - W3x*u; %Right hand side, u_t=-uu_x-u_{xxx}
init = U(x,0.0);
options = odeset('RelTol',2.3e-12,'AbsTol',1e-16);
[t,w] = ode113(@(t,u) RHS_u(t,u),tspan,init,options);
w1 = w(end,:); 
w1 = w1';
norm(w1 - U(x,Tfin),inf)
plot(x,w1,'o',x,U(x,Tfin))