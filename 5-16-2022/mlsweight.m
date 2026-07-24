
function w  = mlsweight (char,x)

global delta;   % size of supports of MLS weight function
global ep; % shape parameter in Gaussian weights.
global weighttype

switch lower(weighttype)
    case ('gauss') % Gaussian weight function
        r = sqrt(sum(x'.^2)')./delta;
        switch lower(char)
            case ('w')
                w=(exp(-ep^2*r.^2)-exp(-ep^2))/(1-exp(-ep^2)); % Gaussian
            case ('wx')
                w=-2*ep^2*x(:,1).*exp(-ep^2*r.^2)/(delta*(1-exp(-ep^2)));
            case ('wy')
                w=-2*ep^2*x(:,2).*exp(-ep^2*r.^2)/(delta*(1-exp(-ep^2)));
            case ('wxx')
                ss = -2*ep^2/(delta*(1-exp(-ep^2)));
                pp = 1-2*ep^2/delta*x(:,1);
                w = ss*pp.*exp(-ep^2*r.^2);
            case ('wyy')
                ss = -2*ep^2/(delta*(1-exp(-ep^2)));
                pp = 1-2*ep^2/delta*x(:,2);
                w = ss*pp.*exp(-ep^2*r.^2);
            case ('wxy')                
                ss = 4*ep^4/(delta^2*(1-exp(-ep^2)));                
                w = ss*x(:,1).*x(:,2).*exp(-ep^2*r.^2);                               
            otherwise
                error ('this derivative of MLS weight is not implemented')
        end
    case ('spline3') % cubic spline weight function
        r = sqrt(sum(x'.^2)');
        s = r/delta;
        n1 = find (s < 0.5 );
        n2= find (s >= 0.5);
        s1 = s(n1); s2 = s(n2);
        switch lower(char)
            case ('w')
                w1 = 2/3 - 4*s1.^2 + 4*s1.^3;
                w2 = 4/3 - 4*s2 + 4*s2.^2 - 4/3*s2.^3;
                w(n1)=w1;
                w(n2)=w2;
            case ('wx')
                w1 = - 8 + 12*s1;
                w2 = - 4./s2 + 8 - 4*s2;
                w(n1)=w1;
                w(n2)=w2;
                w = x(:,1)'/delta^2 .* w;
            case ('wy')
                w1 = - 8 + 12*s1;
                w2 = - 4./s2 + 8 - 4*s2;
                w(n1)=w1;
                w(n2)=w2;
                w = x(:,2)'/delta^2 .* w;
            case ('wxx')
                w1 = - 8 + 12*s1;
                w2 = - 4./s2 + 8 - 4*s2;
                w(n1)=w1;
                w(n2)=w2;
                w = 1/delta^2 .* w;
                n0 = find (abs(s1) <eps);
                w1 = 12./s1;
                w1(n0)=0;
                w2 = 4./s2.^3 - 4./s2;
                ww(n1)=w1;
                ww(n2)=w2;
                w =w +  x(:,1)'.^2/delta^4 .* ww;

            case ('wyy')
                w1 = - 8 + 12*s1;
                w2 = - 4./s2 + 8 - 4*s2;
                w(n1)=w1;
                w(n2)=w2;
                w = 1/delta^2 .* w;
                n0 = find (abs(s1) <eps);
                w1 = 12./s1;
                w1(n0)=0;
                w2 = 4./s2.^3 - 4./s2;
                ww(n1)=w1;
                ww(n2)=w2;
                w =w +  x(:,2)'.^2/delta^4 .* ww;
            case ('wxy')
                w1 = - 8 + 12*s1;
                w2 = - 4./s2 + 8 - 4*s2;
                w(n1)=w1;
                w(n2)=w2;
                w = 1/delta^2 .* w;
                n0 = find (abs(s1) <eps);
                w1 = 12./s1;
                w1(n0)=0;
                w2 = 4./s2.^3 - 4./s2;
                ww(n1)=w1;
                ww(n2)=w2;
                w =w +  x(:,1)'.*x(:,2)'/delta^4 .* ww;

            otherwise
                error ('this derivative of MLS weight is not implemented')
        end
    case ('spline4') % quartic spline weight function
        r = sqrt(sum(x'.^2)');
        s = r/delta;
        switch lower(char)
            case ('w')
                w = 1 - 6*s.^2 + 8*s.^3 - 3*s.^4;
            case ('wx')
                w = x(:,1)/delta^2.* (-12 + 24*s - 12*s.^2);
            case ('wy')
                w = x(:,2)/delta^2.* (-12 + 24*s - 12*s.^2);
            case ('wxx')
                w = x(:,1).^2/delta^4 .* (24./s - 24);
                n0 = find (abs(s) <eps);
                w(n0)=0;
                w = w + 1/delta^2.* (-12 + 24*s - 12*s.^2);
            case ('wyy')
                w = x(:,2).^2/delta^4 .* (24./s - 24);
                n0 = find (abs(s) <eps);
                w(n0)=0;
                w = w + 1/delta^2.* (-12 + 24*s - 12*s.^2);
            case ('wxy')
                w = x(:,1).*x(:,2)/delta^4 .* (24./s - 24);
                n0 = find (abs(s) <eps);
                w(n0)=0;
                
            otherwise
                error ('this derivative of MLS weight is not implemented')
        end
    otherwise
        error ('this type of MLS weight function is not implemented')

end