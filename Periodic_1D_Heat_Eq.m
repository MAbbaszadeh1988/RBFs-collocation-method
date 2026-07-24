%% RBF collocation method for 1D Heat Eq.
clear
clc
format shorte
k = 1;
for dt = [1/10 1/20 1/40 1/80 1/160]
    %% Collocation Points
    n = 20;
    L0 = 0;  L1 = 1;
    h = (L1-L0)/n;
    x = L0:h:L1;
    boundary = find( x==L0 | x==L1 );
    Interior = find( x~=L0 & x~=L1 );
    N = length(x);
    %% Time discrete scheme
    
    T = 1;
    M = floor(T/dt);
    %% Basis Functions and Required Matrices
        Phi = @(r,c) exp(-(c*r).^2);
        drPhi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
        d2rPhi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));
    
    
%     Phi = @(r,c) sqrt( 1 + (c*r).^2 );
%     drPhi = @(r,rx,c) c^2*rx./sqrt( 1 + (c*r).^2 );
%     d2rPhi = @(r,rx,c) (c^2*( 1 + (c*r).^2 ) -c^4*rx.^2)./( ( 1 + (c*r).^2 ).^1.5 );
%     
    rx =  x'-x ;
    D = abs(rx);  % abs(x-x');
    c = 10;
    
    A = Phi(D,c);
    Ax = drPhi(D,rx,c);
    Axx = d2rPhi(D,rx,c);
    %%
    u = @(x,t) exp(-t).*sin(pi.*x);
    f = @(x,t) (-1+pi^2).*exp(-t).*sin(pi.*x);
    
    %%
    A_left  = A - 0.5*dt*Axx;
    A_right = A + 0.5*dt*Axx;
        
    lambda = pinv(A)*u(x',0);
    
    for n = 1:M
        
        F = A_right*lambda + dt*f(x',(n-0.5)*dt);
                
        lambda = pinv(A_left)*F;

        U = A*lambda;

        U(1) = U(end-1); 
        U(end) = U(2);

        lambda = pinv(A)*U;
        
    end
    
    U = A*lambda;
    
    E(k) = norm(U-u(x',T),inf);
    if k~=1
        order(k-1) = log2(E(k-1)/E(k));
    end
    k = k+1;
end
E'
order'