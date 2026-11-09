function R = qstab_eigen_assess(A,reference,eigenvalues,varargin)
%QSTAB_EIGEN_ASSESS Independent reference/multiplicity and eigenpair checks.
% This is a comparison diagnostic, not a replacement eigenvalue algorithm.
% V and D are checked in their ORIGINAL quaternion orientation. Canonical
% representatives are used only for class matching, never in A*V=V*D.
ip=inputParser;
addParameter(ip,'Tolerance',1e-7,@positive);
addParameter(ip,'Info',[],@(x)isempty(x)||(isstruct(x)&&isscalar(x)));
addParameter(ip,'CheckVectors',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Vectors',[]);addParameter(ip,'Diagonal',[]);
addParameter(ip,'Diagonalizable',true,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});o=ip.Results;
[X,xi]=qstab.chi(A);n=xi.quaternionSize(1);
if n~=xi.quaternionSize(2)||numel(reference)~=n
    error('qstab:eigencomparison:Order','A square matrix and n reference classes are required.');
end
R=struct('status','INVALID','ok',false,'reference',reference(:), ...
    'components',[],'canonical',[],'standardForm',false,'spectrum',[], ...
    'solverValid',true,'solverInfo',o.Info,'vectors', ...
    struct('status','NOT_REQUESTED','ok',true),'message','');
try
    [p,z,standard]=representatives(eigenvalues,n,o.Tolerance);
    R.components=p;R.canonical=z;R.standardForm=standard;
    R.spectrum=qstab_eigen_match(z,reference,o.Tolerance);
catch ME
    R.message=ME.message;return
end
if ~isempty(o.Info)
    if ~isfield(o.Info,'isValid')||~isscalar(o.Info.isValid)|| ...
            ~(isnumeric(o.Info.isValid)||islogical(o.Info.isValid))|| ...
            ~isreal(o.Info.isValid)||~isfinite(o.Info.isValid)||o.Info.isValid~=1
        R.solverValid=false;
    end
end
R.status=R.spectrum.status;R.ok=R.spectrum.ok&&R.solverValid;
R.message=R.spectrum.message;
if ~R.solverValid,R.status='UNRESOLVED';R.message='Solver did not mark its result valid.';end
if o.CheckVectors
    R.vectors=checkVectors(X,n,o.Vectors,o.Diagonal,z,o);
    R.ok=R.ok&&R.vectors.ok;
    if ~R.vectors.ok,R.status='EIGENVECTOR_ISSUE';R.message=R.vectors.message;end
end
end
function [p,z,standard]=representatives(E,n,tol)
[a,b,c,d]=qstab.components(E);
if numel(a)~=n||(~isempty(a)&&~isvector(a))
    error('qstab:eigencomparison:Shape','Eigenvalues must be an n-element vector.');
end
p=struct('a',a(:),'b',b(:),'c',c(:),'d',d(:));
r=hypot(hypot(b(:),c(:)),d(:));z=complex(a(:),r);
if any(~isfinite(z)),error('qstab:eigencomparison:Finite','Nonfinite eigenvalue output.');end
s=max(abs(z));if isempty(s)||s==0,s=1;end
standard=all(b(:)/s>=-tol)&&all(abs(c(:))/s<=tol)&&all(abs(d(:))/s<=tol);
end
function R=checkVectors(X,n,V,D,oneOutput,o)
R=struct('status','INVALID','ok',false,'message','', ...
    'V',[],'D',[],'diagonalClasses',[],'diagonalAgreement',[], ...
    'globalResidual',NaN,'columnResiduals',[],'allColumnsNonzero',false, ...
    'complexBasisRank',NaN,'complexOrder',2*n,'fullBasisObserved',false, ...
    'fullBasisRequired',o.Diagonalizable,'rankTolerance',NaN);
try
    [va,vb,vc,vd]=qstab.components(V);[da,db,dc,dd]=qstab.components(D);
    if ~isequal(size(va),[n n])||~isequal(size(da),[n n])
        error('qstab:eigencomparison:VectorShape','V and D must both be n-by-n.');
    end
    R.V=struct('a',va,'b',vb,'c',vc,'d',vd);
    R.D=struct('a',da,'b',db,'c',dc,'d',dd);
    Dparts={da,db,dc,dd};off=0;scaleD=0;
    for k=1:4
        B=Dparts{k};off=max(off,norm(B-diag(diag(B)),'fro'));
        scaleD=max(scaleD,norm(B,'fro'));
    end
    if off>o.Tolerance*max(realmin,scaleD)
        error('qstab:eigencomparison:Diagonal','D is not diagonal.');
    end
    Q=struct('a',diag(da),'b',diag(db),'c',diag(dc),'d',diag(dd));
    [~,z]=representatives(Q,n,o.Tolerance);R.diagonalClasses=z;
    R.diagonalAgreement=qstab_eigen_match(z,oneOutput,o.Tolerance);
    Y=qstab.chi(V);T=qstab.chi(D);
    xs=max([abs(X(:));abs(T(:))]);if isempty(xs)||xs==0,xs=1;end
    ys=max(abs(Y(:)));if isempty(ys),ys=0;end
    if ys==0
        R.message='A zero eigenvector matrix has no valid eigenvectors.';return
    end
    Xn=X/xs;Tn=T/xs;Yn=Y/ys;
    Z=Xn*Yn-Yn*Tn;
    den=norm(Xn,'fro')*norm(Yn,'fro')+norm(Yn*Tn,'fro');
    R.globalResidual=norm(Z,'fro')/max(realmin,den);
    R.columnResiduals=Inf(n,1);nonzero=true;
    for k=1:n
        cols=[k n+k];Yk=Yn(:,cols);vk=norm(Yk,'fro');
        if vk==0,nonzero=false;continue;end
        Yunit=Yk/vk;Tk=Tn(cols,cols);
        R.columnResiduals(k)=norm(Xn*Yunit-Yunit*Tk,'fro')/ ...
            max(realmin,norm(Xn,'fro')+norm(Tk,'fro'));
    end
    R.allColumnsNonzero=nonzero;
    singular=svd(Yn);R.rankTolerance=max(size(Yn))*eps*max(singular);
    R.complexBasisRank=sum(singular>R.rankTolerance);
    R.fullBasisObserved=R.complexBasisRank==2*n;
    R.ok=nonzero&&R.diagonalAgreement.ok&& ...
        isfinite(R.globalResidual)&&R.globalResidual<=o.Tolerance&& ...
        all(isfinite(R.columnResiduals))&&all(R.columnResiduals<=o.Tolerance)&& ...
        (~o.Diagonalizable||R.fullBasisObserved);
    if R.ok
        R.status='CHECKED';R.message='Original oriented eigenpairs and required basis checks passed.';
        if ~o.Diagonalizable,R.message='Eigenpairs checked; a full basis is not required for this defective fixture.';end
    else
        R.status='EIGENVECTOR_ISSUE';
        R.message='A nonzero-column, diagonal, residual, multiplicity or required basis check failed.';
    end
catch ME
    R.message=ME.message;
end
end
function tf=positive(x)
tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0;
end
