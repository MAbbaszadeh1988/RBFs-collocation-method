clear
clc
N = 200;
a = -2;  b = 2;
x = a:(b-a)/N:b;
c=5
hold on
for xi = [-0.5 0 0.5];
    phi = @(x,xi) exp(-(c.*(x-xi)).^2);
    plot(x,phi(x,xi)), drawnow
end