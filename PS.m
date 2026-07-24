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

T = 1;
dt = 1e-4;
Mt = floor(T/dt);


Phi = @(r,c) exp(-(c*r).^2);
drPhi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
d2rPhi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));


D = abs(x-x');
rx =  x'-x ;
A = Phi(D,c);
Ax = drPhi(D,rx,c);
Axx = d2rPhi(D,rx,c);
%%
u = @(x,t) exp(t).*sin(pi.*x);
f = @(x,t) (1+pi^2).*exp(t).*sin(pi.*x);
%%
I = eye(N);

DL = I-0.5*dt*Axx*pinv(A);
DR = I+0.5*dt*Axx*pinv(A);

DL(Boundary,:) = I(Boundary,:);

U = u(x',0);

for n=1:Mt
    t = (n-0.5)*dt;
    F = DR*U + dt*f(x',t);
    F(Boundary) = u(x(Boundary),n*dt);
    U = DL\F;
end
norm(U-u(x',T),inf)

