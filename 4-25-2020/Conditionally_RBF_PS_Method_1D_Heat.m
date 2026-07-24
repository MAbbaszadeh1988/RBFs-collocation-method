clear
clc
format shorte
%%
c = 10;
n = 40;
L0 = 0;  L1 = 1;
h = (L1-L0)/n;
x = L0:h:L1;
Boundary = find(x==L0 | x==L1);
N = length(x);
%%
T = 1;
dt = 1e-4;
Mt = floor(T/dt);
%%
Phi = @(r,c) sqrt( 1 + (c*r).^2 );
drPhi = @(r,rx,c) c^2*rx./sqrt( 1 + (c*r).^2 );
d2rPhi = @(r,rx,c) (c^2*( 1 + (c*r).^2 ) -c^4*rx.^2)./( ( 1 + (c*r).^2 ).^1.5 );

D = abs(x-x');
rx =  x'-x ;
A = Phi(D,c);
Ax = drPhi(D,rx,c);
Axx = d2rPhi(D,rx,c);
%%
m = 10;
P = ones(N,m);
dP = zeros(N,m);
d2P = zeros(N,m);
xp = x';
for j=2:m
    P(:,j) = xp.^(j-1);
    dP(:,j) = (j-1).*xp.^(j-2);
    if j~=2
        d2P(:,j) = (j-1)*(j-2).*xp.^(j-3);
    end
end
%%
u = @(x,t) exp(-t).*(sin(pi.*x)+cos(pi.*x));
f = @(x,t) exp(-t).*(pi^2*cos(pi.*x) + pi^2*sin(pi.*x)) - exp(-t).*(cos(pi.*x) + sin(pi.*x));
%%
M = [A P;P' zeros(m)];
Mx = [Ax dP;P' zeros(m)];
Mxx = [Axx d2P;P' zeros(m)];

Minv = pinv(M);
I = eye(N);
H = Mxx*Minv;

DL = I-0.5*dt*H(1:N,1:N);
DR = I+0.5*dt*H(1:N,1:N);
DL(Boundary,:) = I(Boundary,:);
U = u(x',0);
DLinv = pinv(DL);
for n=1:Mt
    t = (n-0.5)*dt;
    F = DR*U + dt*f(x',t);
    F(Boundary) = u(x(Boundary),n*dt);
    U = DLinv*F;

end
norm(U-u(xp,T),inf)

