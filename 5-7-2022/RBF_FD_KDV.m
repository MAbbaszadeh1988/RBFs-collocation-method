clear all;
close all;
clc;
tic
%%
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
phi    = @(r,c) exp(-(c*r).^2);
drphi  = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
d2rphi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));
d3rphi = @(r,rx,c) 12*c^4.*rx.*exp(-c^2.*r.^2) - 8*c^6*rx.^3.*exp(-c^2*r.^2);

W1x = zeros(M);  % The first-order differential matrix D_{x}
W3x = zeros(M);  % The Third-order differential matrix D_{xxx}
%% Optimal shape parameter S. Sarrah
c = 2;
minK = 1e8;
maxK = 1e10;
dc = 0.1;
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
    B3_local = d3rphi(sqrt((x(i)-x(pn)).^2),x(i)-x(pn),c);
    W1x(i,pn) = pinv(A_local)*B1_local;
    W3x(i,pn) = pinv(A_local)*B3_local;
end
%%
%Zero flux boundary conditions
W1x(1,:) = zeros(size(W1x(1,:))); W1x(end,:) = zeros(size(W1x(1,:)));
W3x(1,:) = zeros(size(W1x(1,:))); W3x(end,:) = zeros(size(W1x(1,:)));
%% Solution of KdV PDE
A = sqrt(1)/sqrt(6); 
L = 1; 
x0 = 0;
U = @(x,t) 3*A^2*sech(A*L*(x - x0/L)/2 - A^3*t/2).^2; %Exact solution
%% 
RHS_u = @(t,u) -u.*(W1x*u) - W3x*u; %Right hand side, u_t=-uu_x-u_{xxx}
init = U(x,0.0);
options = odeset('RelTol',2.3e-12,'AbsTol',1e-16);
[t,w] = ode113(@(t,u) RHS_u(t,u),tspan,init,options);
%%
w1 = w(end,:); 
w1 = w1';
norm(w1 - U(x,Tfin),inf)
plot(x,w1,'o',x,U(x,Tfin))
toc