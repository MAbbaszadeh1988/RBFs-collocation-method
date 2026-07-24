function mat = MLSmat(char,ptrial)
% This function produces MLS shape function and its derivatives
% char = 'a' for shape function itself,
% char = 'ax', 'ay' 1st derivatives in respect to x and y respectively.
% char = 'axx', 'ayy' and 'axy' for 2nd derivatives.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global delta   % size of supports of MLS weight function
global h % fill distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N, dim]=size(ptrial); % N: number of trial points, dim: dimension
mat=zeros(N);
for i=1:N
    X = ones(N,dim)*diag(ptrial(i,:));
    if (dim==1)
        neighbors=delta-abs(X-ptrial)>0;
    else
        neighbors=delta-sqrt(sum(((X-ptrial).^2)')')>0;
    end
    ind=find(neighbors);
    x = X(ind,:)-ptrial(ind,:);
    W= diag(mlsweight ('w', x));
    P=MLSbasisPoly('p',-x/h);
    B=P'*W;
    A=P'*W*P;
    Ainv=A\eye(size(A,1));
    alpha=Ainv*B;
    switch (char)
        case ('a')
            pp= MLSbasisPoly('p',zeros(1,dim));
            mat(i,ind)=pp*alpha;
        case ('ax')
            dw= diag(mlsweight ('wx', x));
            DB=P'*dw;
            DA=DB*P;
            A1=-Ainv*DA;
            DAinv=A1*Ainv;
            pp= MLSbasisPoly('p',zeros(1,dim));
            px= MLSbasisPoly('px',zeros(1,dim));
            ax=px*alpha+pp*(Ainv*DB+DAinv*B);
            mat(i,ind)=ax;
        case ('ay')
            dw= diag(mlsweight ('wy', x));
            DB=P'*dw;
            DA=DB*P;
            A1=-Ainv*DA;
            DAinv=A1*Ainv;
            pp= MLSbasisPoly('p',zeros(1,dim));
            py= MLSbasisPoly('py',zeros(1,dim));
            ay=py*alpha+pp*(Ainv*DB+DAinv*B);
            mat(i,ind)=ay;
        case ('az')
            dw= diag(mlsweight ('wz', x));
            DB=P'*dw;
            DA=DB*P;
            A1=-Ainv*DA;
            DAinv=A1*Ainv;
            pp= MLSbasisPoly('p',zeros(1,dim));
            pz= MLSbasisPoly('pz',zeros(1,dim));
            az=pz*alpha+pp*(Ainv*DB+DAinv*B);
            mat(i,ind)=az;

        case ('axx')
            wx= diag(mlsweight ('wx', x));
            wxx= diag(mlsweight ('wxx', x));
            Bx = P'*wx;
            Bxx=P'*wxx;
            Ax=Bx*P;
            Axx=Bxx*P;
            A1x=Ainv*Ax;
            Ainvx=-A1x*Ainv;
            Ainvxx=2*A1x*A1x*Ainv-Ainv*Axx*Ainv;
            pp= MLSbasisPoly('p',zeros(1,dim));
            px= MLSbasisPoly('px',zeros(1,dim));
            pxx= MLSbasisPoly('pxx',zeros(1,dim));
            Gx=Ainv*Bx+Ainvx*B;
            axx=pxx*alpha+2*px*Gx+pp*(Ainv*Bxx+2*Ainvx*Bx+Ainvxx*B);
            mat(i,ind)=axx;
        case ('ayy')
            wy= diag(mlsweight ('wy', x));
            wyy= diag(mlsweight ('wyy', x));
            By = P'*wy;
            Byy=P'*wyy;
            Ay=By*P;
            Ayy=Byy*P;
            A1y=Ainv*Ay;
            Ainvy=-A1y*Ainv;
            Ainvyy=2*A1y*A1y*Ainv-Ainv*Ayy*Ainv;
            pp= MLSbasisPoly('p',zeros(1,dim));
            py= MLSbasisPoly('py',zeros(1,dim));
            pyy= MLSbasisPoly('pyy',zeros(1,dim));
            Gy=Ainv*By+Ainvy*B;
            ayy=pyy*alpha+2*py*Gy+pp*(Ainv*Byy+2*Ainvy*By+Ainvyy*B);
            mat(i,ind)=ayy;
        case ('azz')
            wx= diag(mlsweight ('wz', x));
            wxx= diag(mlsweight ('wzz', x));
            Bx = P'*wx;
            Bxx=P'*wxx;
            Ax=Bx*P;
            Axx=Bxx*P;
            A1x=Ainv*Ax;
            Ainvx=-A1x*Ainv;
            Ainvxx=2*A1x*A1x*Ainv-Ainv*Axx*Ainv;
            pp= MLSbasisPoly('p',zeros(1,dim));
            px= MLSbasisPoly('pz',zeros(1,dim));
            pxx= MLSbasisPoly('pzz',zeros(1,dim));
            Gx=Ainv*Bx+Ainvx*B;
            axx=pxx*alpha+2*px*Gx+pp*(Ainv*Bxx+2*Ainvx*Bx+Ainvxx*B);
            mat(i,ind)=axx;
        case ('axy')
            wx= diag(mlsweight ('wx', x));
            wy= diag(mlsweight ('wy', x));
            wxy= diag(mlsweight ('wxy', x));
            Bx = P'*wx;
            By = P'*wy;
            Bxy=P'*wxy;
            Ax=Bx*P;
            Ay=By*P;
            Axy=Bxy*P;
            Ainvy=-Ainv*Ay*Ainv;
            Ainvx=-Ainv*Ax*Ainv;
            AYA=Ainv*Ay*Ainv;
            Ainvxy=AYA*Ax*Ainv-Ainv*Axy*Ainv+Ainv*Ax*AYA;
            pp= MLSbasisPoly('p',zeros(1,dim));
            px= MLSbasisPoly('px',zeros(1,dim));
            py= MLSbasisPoly('py',zeros(1,dim));
            pxy= MLSbasisPoly('pxy',zeros(1,dim));
            G1=Ainvy*B+Ainv*By;
            G2=Ainvx*B+Ainv*Bx;
            G3=Ainvxy*B+Ainvx*By+Ainvy*Bx+Ainv*Bxy;
            axy=pxy*alpha+px*G1+ py*G2+pp*G3;
            mat(i,ind)=axy;
        otherwise
            error ('this type of MLS shape function is not implemented')
    end % end switch
end % end for


