clear
clc
%%
f = @(x) exp(x);
a = 0; 
b = 1;
Ng = 10;
%%
[xg,wg]=Legendre_Gauss_Quadrature(Ng,a,b);
% sum = 0;
% for k=1:Ng
%     sum = sum + wg(k)*f(xg(k));
% end
sum = wg'*f(xg);
abs(sum-(exp(1)-1))
abs(quad(f,a,b,1e-14)-(exp(1)-1))
