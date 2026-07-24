% RBF collocation method
clear
clc
format shorte
%% Collocation Points
n = 10;
L0 = 0;  L1 = 1;
h = (L1-L0)/n;
x = L0:h:L1;
boundary = find( x==L0 | x==L1 );
Interior = find( x~=L0 & x~=L1 );
N = length(x);
%% Basis Functions and Required Matrices
Phi = @(r,c) exp(-(c*r).^2);
dxPhi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
dxxPhi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));

% Phi = @(r,c) sqrt( 1 + (c*r).^2 );
% dxPhi = @(r,rx,c) c^2*rx./sqrt( 1 + (c*r).^2 );
% dxxPhi = @(r,rx,c) (c^2*( 1 + (c*r).^2 ) -c^4*rx.^2)./( ( 1 + (c*r).^2 ).^1.5 );


rx =  x'-x ;
D = abs(rx);  % abs(x-x');
c = 20;

Int = Phi(D,c);
D1 = dxPhi(D,rx,c);
D2 = dxxPhi(D,rx,c);
%%
ue = @(x) sin(pi.*x);
f = @(x) -pi^2.*sin(pi.*x)+sin(ue(x));
F = @(u) sin(u);
Fp = @(u) cos(u);
%%

u = x'.^2+1;
u(1) = ue(L0);
u(end) = ue(L1);

count  = 1;

while count==1
    
    H = D2 + Fp(u).*Int;
    H(boundary,:) = Int(boundary,:);
    
    Fr = -D2*u - F(u) + f(x');
    Fr(boundary) = ue(boundary);
    
    Lambda = pinv(H)*Fr;
    
    deltau = Int*Lambda;
    
    u = u + deltau;
    
    if norm(deltau)<1e-12
        count = 0;
    end
    norm(u-ue(x'),inf)
end

norm(u-ue(x'),inf)
plot(x,u,x,ue(x))