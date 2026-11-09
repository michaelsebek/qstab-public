function [lambda,info] = eig(A,varargin)
%EIG Canonical representatives of right quaternion eigenvalue classes.
%   LAMBDA = qstab.eig(A) returns n complex representatives with nonnegative
%   imaginary parts for an n-by-n quaternion matrix. Repeated entries are
%   retained: the vector has algebraic multiplicity n, not just distinct
%   classes. A real eigenvalue of the complex adjoint has doubled
%   multiplicity, and its copies are paired explicitly.
%   [LAMBDA,INFO] reports the full 2n spectrum, conjugate pairing, numerical
%   residuals and diagnostic class groups. Check INFO.isValid before using
%   the representatives. Grouping by ClassTolerance is numerical, not an
%   exact assertion of multiple eigenvalues.
%   Options: 'PairingTolerance' (default 1e-7), 'ClassTolerance' (1e-7).
%   'Tolerance',T sets both unless overridden by a subsequent option.
%
%   Examples (ascending coefficients):
%   [lambda,I] = qstab.eig(diag([-1+1i,-1+2*1i,-1+3*1i]));
%   disp(lambda); disp(I.isValid);
%
%   Interpretation:
%   This is the RIGHT spectrum A*v=v*lambda. The second output is diagnostics, NOT eigenvectors. Use eigpairs for V,D. Do not replace the existing qeig function or its two-output convention.
%
%   Reference: docs/reference/eig.html.
%   See also qstab.about.

[pairTol,classTol]=localOptions(varargin{:});
[X,xi]=qstab.chi(A);
if xi.quaternionSize(1)~=xi.quaternionSize(2)
    error('qstab:NotSquareMatrix','Right eigenvalues require a square quaternion matrix.');
end
[V,D]=builtin('eig',X);
full=diag(D);
[lambda,pairs,pairResidual,valid,message]=localPair(full,pairTol);
if isempty(X)
    residual=0;
else
    VD=V*D;
    residual=norm(X*V-VD,'fro')/max(realmin,norm(X,'fro')*norm(V,'fro')+norm(VD,'fro'));
end
[groupReps,groupMult,groupIndex]=localGroups(lambda,classTol);
isValid=valid && ~isnan(residual) && ~isinf(residual) && residual<=max(pairTol,1e3*eps);
if valid && ~isValid
    message='The embedded eigenvector residual is too large or nonfinite.';
end
info=struct('status','ok','isValid',isValid,'message',message, ...
    'fullSpectrum',full,'pairedIndices',pairs, ...
    'pairingResidual',pairResidual,'pairingTolerance',pairTol, ...
    'pairingValid',valid,'embeddingResidual',xi.embeddingResidual, ...
    'eigenResidual',residual,'classRepresentatives',groupReps, ...
    'classMultiplicities',groupMult,'classIndex',groupIndex, ...
    'classTolerance',classTol,'algebraicCount',numel(lambda), ...
    'quaternionOrder',xi.quaternionSize(1),'spectrumComputed',true, ...
    'assurance','floating-point');
if ~isValid, info.status='undetermined'; end
end

function [pairTol,classTol]=localOptions(varargin)
pairTol=1e-7; classTol=1e-7;
if mod(numel(varargin),2)~=0
    error('qstab:InvalidOptions','Options must be name/value pairs.');
end
for k=1:2:numel(varargin)
    name=lower(char(varargin{k}));
    switch name
        case 'tolerance', pairTol=varargin{k+1}; classTol=varargin{k+1};
        case 'pairingtolerance', pairTol=varargin{k+1};
        case 'classtolerance', classTol=varargin{k+1};
        otherwise, error('qstab:UnknownOption','Unknown option: %s.',name);
    end
end
if ~localTolerance(pairTol) || ~localTolerance(classTol)
    error('qstab:InvalidTolerance','Tolerances must be finite positive real scalars.');
end
end

function tf=localTolerance(x)
tf=isnumeric(x) && isreal(x) && isscalar(x) && ~isnan(x) && ~isinf(x) && x>0;
end

function [lambda,pairs,residual,valid,message]=localPair(full,tol)
full=full(:); count=numel(full); n=count/2;
if mod(count,2)~=0
    error('qstab:InvalidEmbeddingSize','An adjoint spectrum must have even size.');
end
lambda=nan(n,1); pairs=nan(n,2); residual=inf(n,1);
infIdx=find(isinf(full)); finiteIdx=find(~isinf(full) & ~isnan(full));
if any(isnan(full) & ~isinf(full)) || mod(numel(infIdx),2) || mod(numel(finiteIdx),2)
    valid=false; message='The full spectrum cannot be partitioned into finite conjugate pairs and even infinite pairs.'; return
end
row=0;
while ~isempty(finiteIdx)
    [~,p]=max(abs(imag(full(finiteIdx))));
    i=finiteIdx(p); finiteIdx(p)=[];
    distances=abs(full(finiteIdx)-conj(full(i)))./max(1,max(abs(full(finiteIdx)),abs(full(i))));
    [r,p]=min(distances); j=finiteIdx(p); finiteIdx(p)=[];
    row=row+1;
    z=(full(i)+conj(full(j)))/2;
    lambda(row)=complex(real(z),abs(imag(z)));
    pairs(row,:)=[i,j]; residual(row)=r;
end
for k=1:2:numel(infIdx)
    row=row+1; lambda(row)=Inf; pairs(row,:)=infIdx(k:k+1).'; residual(row)=0;
end
if ~isempty(lambda)
    [~,order]=sortrows([real(lambda),imag(lambda)],[1,2]);
    lambda=lambda(order); pairs=pairs(order,:); residual=residual(order);
end
valid=all(residual<=tol);
if valid
    message='Conjugate pairs retained with algebraic multiplicity.';
else
    message='Conjugate pairing exceeds tolerance; the canonical representatives are unresolved.';
end
end

function [reps,mult,index]=localGroups(lambda,tol)
reps=zeros(0,1); mult=zeros(0,1); index=zeros(numel(lambda),1);
for k=1:numel(lambda)
    z=lambda(k); match=[];
    if isinf(z)
        match=find(isinf(reps),1);
    elseif ~isnan(z)
        match=find(abs(reps-z)<=tol*max(1,max(abs(reps),abs(z))),1);
    end
    if isempty(match)
        reps(end+1,1)=z; mult(end+1,1)=0; match=numel(reps);
    end
    mult(match)=mult(match)+1; index(k)=match;
end
end
