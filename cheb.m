function [D,x] = cheb(N,a,b)
if N==0, D=0; x=1; return, end
x = ((b-a)/2).*cos(pi*(0:N)./N)' + ((b+a)/2);
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D = (c*(1./c)')./(dX+(eye(N+1))); % off-diagonal entries
D = D - diag(sum(D')); % diagonal entries
end