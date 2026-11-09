function [L, info] = lyapoperator(A)
%LYAPOPERATOR Real matrix of X -> A*X + X*A' on quaternion Hermitian X.
%   L = qstab.lyapoperator(A) has size (2*n^2-n)-by-(2*n^2-n).
%   A is a square native quaternion/numeric array, or a component struct
%   with real fields a,b,c,d. This is NOT a quaternion Kronecker sum.
%
%   Coordinates use a Frobenius-orthonormal quaternion Hermitian basis:
%   real diagonal entries first; then sqrt(2)*[a,b,c,d] of each upper
%   triangular entry, in row order. INFO documents this coordinate order.
%   The construction uses the complex adjoint only for multiplication.
%
%   Examples (ascending coefficients):
%   [L,I] = qstab.lyapoperator(-eye(2));
%   disp(L);
%
%   Interpretation:
%   The output is a REAL operator on quaternion Hermitian matrices, not a quaternionic Kronecker sum. Its eigenvalue interpretation is part of the directional test; do not substitute a naive noncommutative Kronecker identity.
%
%   Reference: docs/reference/lyapoperator.html.
%   See also qstab.about.

    [a,b,c,d] = qstab.components(A);
    if ndims(a) > 2 || isempty(a) || size(a,1) ~= size(a,2)
        error('qstab:lyapoperator:Shape', 'A must be a nonempty square matrix.');
    end
    n = size(a,1);
    H = qstab.chi(A);
    dim = 2*n*n-n;
    L = zeros(dim,dim);
    pairs = zeros(n*(n-1)/2,2);
    k = 0;
    for i=1:n
        for j=i+1:n
            k=k+1;
            pairs(k,:)=[i,j];
        end
    end
    defect=0;
    for col=1:dim
        E = zeros(n,n,4);
        if col<=n
            E(col,col,1)=1;
        else
            pair=floor((col-n-1)/4)+1;
            part=mod(col-n-1,4)+1;
            i=pairs(pair,1); j=pairs(pair,2);
            E(i,j,part)=1/sqrt(2);
            if part==1
                E(j,i,part)=1/sqrt(2);
            else
                E(j,i,part)=-1/sqrt(2);
            end
        end
        Z=E(:,:,1)+1i*E(:,:,2); W=E(:,:,3)+1i*E(:,:,4);
        Ec=[Z,W;-conj(W),conj(Z)];
        Y=H*Ec+Ec*H';
        Y=(Y+Y')/2;
        Z=(Y(1:n,1:n)+conj(Y(n+1:end,n+1:end)))/2;
        W=(Y(1:n,n+1:end)-conj(Y(n+1:end,1:n)))/2;
        Yproj=[Z,W;-conj(W),conj(Z)];
        defect=max(defect,norm(Y-Yproj,'fro'));
        vec=zeros(dim,1);
        vec(1:n)=real(diag(Z));
        for pair=1:size(pairs,1)
            i=pairs(pair,1); j=pairs(pair,2);
            vec(n+(pair-1)*4+(1:4))=sqrt(2)*...
                [real(Z(i,j));imag(Z(i,j));real(W(i,j));imag(W(i,j))];
        end
        L(:,col)=vec;
    end
    info=struct('matrixOrder',n,'dimension',dim,'pairs',pairs,...
        'coordinateOrder','diagonal, then sqrt(2)*(a,b,c,d) of upper entries',...
        'operator','A*X + X*A''','structureResidual',defect,...
        'assurance','floating-point');
end
