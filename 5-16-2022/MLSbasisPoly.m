function P = MLSbasisPoly(char,xx)
global h; % fill distance
global m;  % polynomaial degree

[n,dim]=size(xx); % dim = dimension of proplem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dim =1
if dim==1
    x =xx;
    switch (char)
        case ('p')
            P = ones(n,1);
            for i = 1:m
                P = [P x.^i];
            end
        case('px')

            P = zeros(n,1);
            for i = 1:m
                P = [P i*x.^(i-1)/h];
            end

        case('pxx')

            P = [zeros(n,1) zeros(n,1)];
            for i = 2:m
                P = [P i*(i-1)*x.^(i-2)/h^2];
            end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dim =2
elseif dim ==2
    x = xx(:,1); y = xx(:,2);
    z0 = zeros(n,1); o1 = ones(n,1);
    switch (char)
        case ('p')
            if (m==1)
                P = [o1 x y];
            elseif (m==2)
                P = [o1 x y x.^2 x.*y y.^2];
            elseif (m==3)
                P = [o1 x y x.^2 x.*y y.^2 x.^3 x.^2.*y x.*y.^2 y.^3];
            elseif (m==4)
                P = [o1 x y x.^2 x.*y y.^2 x.^3 x.^2.*y x.*y.^2 y.^3 ...
                    x.^4 x.^3.*y x.^2.*y.^2 x.*y.^3 y.^4];
            else
                error('Try for lower degrees for basis polynomial')
            end

        case ('px') % derivative with respect to x
            if (m==1)
                P = 1/h*[z0 o1 z0];
            elseif (m==2)
                P = 1/h*[z0 o1 z0 2*x y z0];
            elseif (m==3)
                P = 1/h*[z0 o1 z0 2*x y z0 3*x.^2 2*x.*y y.^2 z0];
            elseif (m==4)
                P = 1/h*[z0 o1 z0 2*x y z0 3*x.^2 2*x.*y y.^2 z0 ...
                    4*x.^3 3*x.^2.*y 2*x.^2.*y.^2 y.^3 z0];
            else
                error('Try for lower degrees for basis polynomial')
            end

        case ('py') % derivative with respect to y
            if (m==1)
                P = 1/h*[z0 z0 o1];
            elseif (m==2)
                P = 1/h*[z0 z0 o1 z0 x 2*y];
            elseif (m==3)
                P = 1/h*[z0 z0 o1 z0 x 2*y z0 x.^2 2*x.*y 3*y.^2];
            elseif (m==4)
                P = 1/h*[z0 z0 o1 z0 x 2*y z0 x.^2 2*x.*y 3*y.^2 ...
                    z0 x.^3 2*x.^2.*y 3*x.*y.^2 4*y.^3];
            else
                error('Try for lower degrees for basis polynomial')
            end

        case ('pxx') % 2d derivative with respect to x
            if (m==1)
                P =1/h^2*[z0 z0 z0];
            elseif (m==2)
                P =1/h^2*[z0 z0 z0 2*o1 z0 z0];
            elseif (m==3)
                P =1/h^2*[z0 z0 z0 2*o1 z0 z0 6*x 2*y z0 z0];
            elseif (m==4)
                P =1/h^2*[z0 z0 z0 2*o1 z0 z0 6*x 2*y z0 z0 12*x.^2 6*x.*y 2*y.^2 z0 z0];
            else
                error('Try for lower degrees for basis polynomial')
            end

        case ('pyy') % 2nd derivative with respect to y

            if (m==1)
                P =1/h^2*[z0 z0 z0];
            elseif (m==2)
                P =1/h^2*[z0 z0 z0 z0 z0 2*o1];
            elseif (m==3)
                P =1/h^2*[z0 z0 z0 z0 z0 2*o1 z0 z0 2*x 6*y];
            elseif (m==4)
                P =1/h^2*[z0 z0 z0 z0 z0 2*o1 z0 z0 2*x 6*y z0 z0 2*x.^2 6*x.*y 12*y.^2];
            else
                error('Try for lower degrees for basis polynomial')
            end

        case ('pxy') % 2nd derivative with respect to y

            if (m==1)
                P =1/h^2*[z0 z0 z0];
            elseif (m==2)
                P = 1/h^2*[z0 z0 z0 z0 o1 z0];
            elseif (m==3)
                P = 1/h^2*[z0 z0 z0 z0 o1 z0 z0 2*x 2*y z0];
            elseif (m==4)
                P = 1/h^2*[z0 z0 z0 z0 o1 z0 z0 2*x 2*y z0 ...
                    z0 3*x.^2 4*x.^2.*y 3*y.^2 z0];
            else
                error('Try for lower degrees for basis polynomial')
            end

        otherwise
            error('Basis type not implemented')
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dim =3
elseif (dim == 3)
    x = xx(:,1); y = xx(:,2); z = xx(:,3);
    z0 = zeros(n,1); o1 = ones(n,1);
    switch (char)
        case ('p')
            if (m==1)
                P = [o1 x y z];
            elseif (m==2)
                P = [o1 x y z x.^2 y.^2 z.^2 x.*y x.*z y.*z];
            elseif (m==3)
                P = [o1 x y z x.^2 y.^2 z.^2 x.*y x.*z y.*z x.^3 y.^3 z.^3 ...
                    x.^2.*y x.^2.*z x.*y.^2 x.*z.^2 y.^2.*z y.*z.^2 x.*y.*z];
            else
                error('Try for lower degrees for basis polynomial')
            end
        otherwise
            error('Basis type not implemented')
    end
else
    error('Polynomials in dimensions more than 2 are not implemented yet')
end


