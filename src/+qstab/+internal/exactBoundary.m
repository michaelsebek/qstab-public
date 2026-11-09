function W = exactBoundary(kind,X,domain)
%EXACTBOUNDARY Sufficient structural witnesses for strict stability failure.
% A tolerance-band root or a rounded zero in a norm polynomial is NEVER
% used. Only signed permutations, comparisons and exact cancellation pairs
% on finite input components enter these decisions. Delays have already
% been reversed at their declared grade. This is not a complete boundary
% solver. Missing evidence deliberately leaves an unresolved result.
W=[]; ct=strcmp(domain,'continuous');
if strcmp(kind,'matrix')
    [a,b,c,d]=qstab.components(X); n=size(a,1);
    nz = a~=0 | b~=0 | c~=0 | d~=0;
    triangular = ~any(any(tril(nz,-1))) || ~any(any(triu(nz,1)));
    % For a triangular quaternion matrix its right spectrum is the union
    % of the diagonal similarity classes. An isolated column or row (only
    % a diagonal entry) also gives a scalar diagonal block by permutation.
    for j=1:n
        if ~triangular && (nnz(nz(:,j))>nz(j,j) && nnz(nz(j,:))>nz(j,j)),continue;end
        q=double([a(j,j);b(j,j);c(j,j);d(j,j)]);
        exact = (ct && q(1)==0) || (~ct && nnz(q)==1 && max(abs(q))==1);
        if exact
            W=localRecord(domain,'diagonalBlock',q,j, ...
                'An exact diagonal-block right eigenvalue lies on the strict stability boundary.');
            return
        end
    end
    return
end
if strcmp(kind,'polynomialmatrix')
    % The original quaternion diagonal polynomial of a common triangular
    % family of coefficient matrices is a factor of the embedded pencil
    % determinant. The same holds for a row/column-isolated scalar block.
    n=X.n; m=2*n; g=X.grade; nz=false(n); diagC=zeros(4,g+1,n);
    for t=0:g
        if t==g, Q=X.B(1:m,1:m);
        else, Q=-X.A(1:m,(g-1-t)*m+(1:m));end
        Z=Q(1:n,1:n); J=Q(1:n,n+1:m);
        nz=nz | Z~=0 | J~=0;
        diagC(:,t+1,:)=reshape([real(diag(Z)).';imag(diag(Z)).'; ...
            real(diag(J)).';imag(diag(J)).'],4,1,n);
    end
    triangular=~any(any(tril(nz,-1))) || ~any(any(triu(nz,1)));
    for j=1:n
        if ~triangular && nnz(nz(:,j))>nz(j,j) && nnz(nz(j,:))>nz(j,j),continue;end
        W=qstab.internal.exactBoundary('polynomial',diagC(:,:,j),domain);
        if ~isempty(W)
            W.diagonalIndex=j; W.reason=['Exact diagonal polynomial block: ',W.reason];return
        end
    end
    return
end
C=X; g=size(C,2)-1;
if g<1 || all(C(:,end)==0),return;end
if ct && all(C(:,1)==0)
    W=localRecord(domain,'zeroConstantTerm',zeros(4,1),[], ...
        'The nonzero polynomial has an exactly zero constant coefficient and therefore a zero at s=0.');
    return
end
% Test actual right evaluation at +1,-1 (DT) and the six imaginary units
% (CT or DT). Right multiplication by a basis unit is a signed permutation.
powers=0:g;
if ~ct
    T=C;
    for sign=[1 -1]
        if sign<0,T(:,2:2:end)=-T(:,2:2:end);end
        if localCancels(T)
            W=localRecord(domain,'signedCoefficientCancellation',[sign;0;0;0],[], ...
                'Exact cancellation at a real unit root confirms failure of strict Schur stability.');
            return
        end
    end
end
for u=1:3
    switch u
        case 1, Cu=[-C(2,:); C(1,:); C(4,:);-C(3,:)];
        case 2, Cu=[-C(3,:);-C(4,:); C(1,:); C(2,:)];
        case 3, Cu=[-C(4,:); C(3,:);-C(2,:); C(1,:)];
    end
    T=C; T(:,2:2:end)=Cu(:,2:2:end); negate=mod(powers,4)>=2;
    T(:,negate)=-T(:,negate);
    for sign=[1 -1]
        if sign<0,T(:,2:2:end)=-T(:,2:2:end);end
        if localCancels(T)
            q=zeros(4,1);q(u+1)=sign;
            W=localRecord(domain,'basisUnitCancellation',q,[], ...
                'Exact cancellation at an imaginary basis unit exhibits a boundary zero.');
            return
        end
    end
end
if g==1
    % p=q0+q1*v. Re(-q1^-1*q0)=0 if the real coefficient vectors
    % have disjoint support. In DT equality of the multisets of component
    % magnitudes implies exactly equal norms, even at extreme scales.
    if ct && ~any(C(:,1)~=0 & C(:,2)~=0)
        W=localRecord(domain,'orthogonalLinearCoefficients',[],[], ...
            'Disjoint coefficient-component support gives an exactly imaginary linear root.');
    elseif ~ct && isequal(sort(abs(C(:,1))),sort(abs(C(:,2))))
        W=localRecord(domain,'equalLinearCoefficientNorms',[],[], ...
            'Matching component magnitudes prove exactly equal coefficient norms and a unit-modulus linear root.');
    end
end
end

function tf=localCancels(T)
tf=true;
for nu=1:4
    x=T(nu,:); p=x(x>0); n=-x(x<0);
    if ~isequal(sort(p(:)),sort(n(:))),tf=false;return;end
end
end

function W=localRecord(domain,kind,point,index,reason)
W=struct('domain',domain,'kind',kind,'quaternionComponents',point, ...
    'diagonalIndex',index,'variableMeaning','forward test variable', ...
    'exactForSuppliedComponents',true,'reason',reason);
end
