clear all
clc
clear all
format shorte
global m % degree of polynomial basis functions (order of approximation)
global h delta % h: fill distance of trial points, delta: support size of MLS weight
global L; % Omega = [0,L]x[0,L]
global weighttype;   % type of MLS weight function
global ep % shape parameter in Gaussian weights
%%
weighttype='gauss'; % weight function, here 'gauss', 'spline3' and 'spline4'
m=4;  % degree of polymomial basis function
L=1; %Omega = [0,L]^2
%%
n=20;
h=L/n;
%%
T = 1;  dt = 1e-3;  Nt = T/dt;
%%
[xt,yt]=meshgrid(-L:h:L,-L:h:L);
ptrial=[xt(:) yt(:)];
[boundary] = find(xt(:)==-L | xt(:)==L | yt(:)==-L | yt(:)==L);
delta=2*m*h; % the size of MLS supports
ep = 10; % shape parameter in Gaussian weight (experimentally)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Imat = MLSmat('axx',ptrial) + MLSmat('ayy',ptrial);

Amat = MLSmat('a',ptrial);

I = eye(size(Amat));

uexact =@(x,y,t) t^2.*cos(pi.*x).*cos(pi.*y);
f = @(x,y,t) (2*t+2*t^2*pi^2).*cos(pi.*x).*cos(pi.*y);

A_left  = I - (dt/2).*Imat;
A_right = I + (dt/2).*Imat;

A_left(boundary,:) = I(boundary,:);

U = uexact(xt(:),yt(:),0);

for k=1:Nt


    % Calculation
%     workbar(k*dt,'Performing Calclations...','Progress')

    t = (k-0.5)*dt
    F = A_right*U + dt.*f(xt(:),yt(:),t);
    F(boundary) = uexact(xt(boundary),yt(boundary),k*dt);
    U = A_left \F;
end
clc
norm(U-uexact(xt(:),yt(:),T),inf)
figure, pcolor(xt,yt,reshape(U,size(xt))), shading interp
figure, pcolor(xt,yt,reshape(abs(U-uexact(xt(:),yt(:),T)),size(xt))), shading interp