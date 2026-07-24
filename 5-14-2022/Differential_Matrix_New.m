function [Phi]=Differential_Matrix_New(gpos,x,y,c,N,mbasis,char)

rx = x - x';
ry = y - y';

r = sqrt( rx.^2 + ry.^2);

A   = sqrt(1+c^2.*r.^2);  
%%
P0t = zeros(mbasis,N);
P0t(1,:) = 1;
P0t(2,:) = x;
P0t(3,:) = y;
P0t(4,:) = x.^2;
P0t(5,:) = x.*y;
P0t(6,:) = y.^2;

M = [A P0t';P0t zeros(mbasis)];

Minv = pinv(M);
%%
R = sqrt(rx.^2+ry.^2);

Phi = zeros(N,N+mbasis);

switch  (char)
    case ('Int')
        R0 = zeros(N,N+mbasis);
        R0(1:N,1:N) = sqrt(1+c^2.*R.^2);
            R0(1:N,N+1) = 1;
            R0(1:N,N+2) = gpos(1,1:N);
            R0(1:N,N+3) = gpos(2,1:N);
            R0(1:N,N+4) = gpos(1,1:N).^2;
            R0(1:N,N+5) = gpos(1,1:N).*gpos(2,1:N);
            R0(1:N,N+6) = gpos(2,1:N).^2;
        Phi = R0*Minv;
    case ('ax')
        R0x = zeros(N,N+mbasis);
        R0x(1:N,1:N) = c^2.*(rx)./sqrt(1+c^2.*R.^2);
        for i=1:N
            R0x(i,N+1) = 0;
            R0x(i,N+2) = 1;
            R0x(i,N+3) = 0;
            R0x(i,N+4) = 2*gpos(1,i);
            R0x(i,N+5) = gpos(2,i);
            R0x(i,N+6) = 0;
        end
        Phi = R0x*Minv;
    case ('ay')
        R0y = zeros(N,N+mbasis);
        R0y(1:N,1:N) = c^2.*(ry)./sqrt(1+c^2.*R.^2);
        for i=1:N
            R0y(i,N+1) = 0;
            R0y(i,N+2) = 0;
            R0y(i,N+3) = 1;
            R0y(i,N+4) = 0;
            R0y(i,N+5) = gpos(1,i);
            R0y(i,N+6) = 2*gpos(2,i);
        end
        Phi = R0y*Minv;
    case ('axx')
        R0xx = zeros(N,N+mbasis);
        R0xx(1:N,1:N) = (c^2.*(1+c^2.*R.^2)-c^4.*(rx).^2)./(((1+c^2.*R.^2)).^1.5);
        for i=1:N
            R0xx(i,N+1) = 0;
            R0xx(i,N+2) = 0;
            R0xx(i,N+3) = 0;
            R0xx(i,N+4) = 2;
            R0xx(i,N+5) = 0;
            R0xx(i,N+6) = 0;
        end
        Phi = R0xx*Minv;
    case ('ayy')
        R0yy = zeros(N,N+mbasis);
        R0yy(1:N,1:N) = (c^2.*(1+c^2.*R.^2)-c^4.*(ry).^2)./(((1+c^2.*R.^2)).^1.5);
        for i=1:N
            R0yy(i,N+1) = 0;
            R0yy(i,N+2) = 0;
            R0yy(i,N+3) = 0;
            R0yy(i,N+4) = 0;
            R0yy(i,N+5) = 0;
            R0yy(i,N+6) = 2;
        end
        Phi = R0yy*Minv;
end
