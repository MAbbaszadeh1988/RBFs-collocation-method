clear
clc
n = 80;
L0 = -2;
L1 = 2;
h = (L1-L0)/n;
x = L0:h:L1;
[X,Y] = meshgrid(x);

f = @(x,y) sin(x+y);

% plot3(X(:),Y(:),f(X(:),Y(:)),'o')

% mesh(X,Y,f(X,Y))

surf(X,Y,f(X,Y)), shading interp

% pcolor(X,Y,f(X,Y)), shading interp
