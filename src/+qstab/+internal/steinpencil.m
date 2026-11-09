function [L0,L1,L2,info] = steinpencil(A0,A1,E0,E1)
%STEINPENCIL Real quadratic generalized Stein operator on Hermitian H^n.
% E(beta)*X*E(beta)' - A(beta)*X*A(beta)', E=E0+beta*E1, A=A0+beta*A1.
% Inputs are structure-preserving complex adjoints (2n)-by-(2n).
% Coordinates match qstab.lyapoperator: real diagonal then sqrt(2)*a,b,c,d.
n=size(A0,1)/2;dim=2*n*n-n;L0=zeros(dim);L1=L0;L2=L0;
pairs=zeros(n*(n-1)/2,2);k=0;
for i=1:n,for j=i+1:n,k=k+1;pairs(k,:)=[i,j];end;end
structureResidual=0;
for col=1:dim
    Q=zeros(n,n,4);
    if col<=n,Q(col,col,1)=1;
    else
        pair=floor((col-n-1)/4)+1;part=mod(col-n-1,4)+1;
        i=pairs(pair,1);j=pairs(pair,2);Q(i,j,part)=1/sqrt(2);
        if part==1,Q(j,i,part)=1/sqrt(2);else,Q(j,i,part)=-1/sqrt(2);end
    end
    Z=Q(:,:,1)+1i*Q(:,:,2);W=Q(:,:,3)+1i*Q(:,:,4);X=[Z,W;-conj(W),conj(Z)];
    Y0=E0*X*E0'-A0*X*A0';
    Y1=E0*X*E1'+E1*X*E0'-A0*X*A1'-A1*X*A0';
    Y2=E1*X*E1'-A1*X*A1';
    [L0(:,col),r0]=coordinates(Y0,n,pairs,dim);
    [L1(:,col),r1]=coordinates(Y1,n,pairs,dim);
    [L2(:,col),r2]=coordinates(Y2,n,pairs,dim);
    structureResidual=max([structureResidual,r0,r1,r2]);
end
info=struct('matrixOrder',n,'dimension',dim,'pairs',pairs, ...
    'operator','E(beta)*X*E(beta)'' - A(beta)*X*A(beta)''', ...
    'coordinateOrder','diagonal, then sqrt(2)*(a,b,c,d) of upper entries', ...
    'structureResidual',structureResidual,'validated',false);
end
function [v,defect]=coordinates(Y,n,pairs,dim)
Y=(Y+Y')/2;
Z=(Y(1:n,1:n)+conj(Y(n+1:end,n+1:end)))/2;
W=(Y(1:n,n+1:end)-conj(Y(n+1:end,1:n)))/2;
projected=[Z,W;-conj(W),conj(Z)];defect=norm(Y-projected,'fro');
v=zeros(dim,1);v(1:n)=real(diag(Z));
for k=1:size(pairs,1)
    i=pairs(k,1);j=pairs(k,2);
    v(n+4*(k-1)+(1:4))=sqrt(2)*[real(Z(i,j));imag(Z(i,j));real(W(i,j));imag(W(i,j))];
end
end
