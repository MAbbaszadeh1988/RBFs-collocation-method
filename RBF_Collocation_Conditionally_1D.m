clear
clc
for i = 1:50
    c =   sqrt(i);
    n = 1000;
    L0 = 0;  L1 = 1;
    h = (L1-L0)/n;
    x = L0:h:L1;
    N = length(x);
    D = abs(x-x');
    %%
    Phi = @(r,c) sqrt( 1 + (c*r).^2 );
    drPhi = @(r,rx,c) c^2*rx./sqrt( 1 + (c*r).^2 );
    d2rPhi = @(r,rx,c) (c^2*( 1 + (c*r).^2 ) -c^4*rx.^2)./( ( 1 + (c*r).^2 ).^1.5 );
    
    rx =  x'-x ;
    A = Phi(D,c);
    Ax = drPhi(D,rx,c);
    Axx = d2rPhi(D,rx,c);
    
    %%
    m = 2;
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
    u = @(x) sin(pi.*x);
    f = @(x) (1-pi^2).*sin(pi.*x)+pi.*cos(pi.*x);
    %%
    M = [A P;P' zeros(m)];
    Mx = [Ax dP;zeros(m,N) zeros(m)];
    Mxx = [Axx d2P;zeros(m,N) zeros(m)];
    
%     H = [Axx+Ax+A d2P+dP+P ; P' zeros(m)];
    
    H =  Mxx + Mx + M;
    
    H(1,:) = M(1,:);
    H(N,:) = M(N,:);
    
    F = [f(xp);zeros(m,1)];
    F(1) = u(L0);
    F(N) = u(L1);
    lambda = H\F;
    U = M*lambda;
    %%
    E(i) = norm(U(1:N)-u(x'),inf);
    Shape(i) = c;
end
loglog(Shape,E,'-.')