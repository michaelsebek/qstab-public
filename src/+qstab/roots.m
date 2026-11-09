function [lambda,info] = roots(P,varargin)
%ROOTS Algebraic zero classes of a quaternion polynomial or matrix polynomial.
%   [LAMBDA,INFO] = qstab.roots(P) uses an embedded generalized companion
%   pencil. Ascending coefficients are required. Scalar coefficient vectors,
%   square coefficient tensors, coefficient cells and a,b,c,d structs are
%   accepted. LAMBDA is a complex column with nonnegative imaginary parts;
%   multiplicities are retained, including real and infinite classes.
%   This computes algebraic classes, not orientation-specific evaluation
%   zeros. INFO.fullSpectrum contains the complete embedded spectrum.
%   INFO.isValid reports successful conjugate pairing and a numerical
%   regularity witness. Singular polynomial matrices are not solved here.
%   Options: 'Grade',G; 'Variable','s'|'z'|'zinv' (metadata only);
%   'PairingTolerance',T (1e-7); 'ClassTolerance',T (1e-7); 'Tolerance',T
%   sets both tolerances. Roots remain in the input variable, also for zinv.
%
%   Examples (ascending coefficients):
%   [lambda,I] = qstab.roots([2 3 1],'Variable','s');
%   disp(lambda); disp(I.isValid);
%
%   Interpretation:
%   These are canonical representatives of algebraic zero classes, with multiplicity. They are not orientation-specific evaluation zeros, not the old qroots contract, and not left eigenvalues. A delay output stays in its input variable.
%
%   Reference: docs/reference/roots.html.
%   See also qstab.about.

grade=[]; variable='s'; pairTol=1e-7; classTol=1e-7;
if mod(numel(varargin),2)~=0
    error('qstab:InvalidOptions','Options must be name/value pairs.');
end
for k=1:2:numel(varargin)
    name=lower(char(varargin{k}));
    switch name
        case 'grade', grade=varargin{k+1};
        case 'variable', variable=lower(char(varargin{k+1}));
        case 'tolerance', pairTol=varargin{k+1}; classTol=varargin{k+1};
        case 'pairingtolerance', pairTol=varargin{k+1};
        case 'classtolerance', classTol=varargin{k+1};
        otherwise, error('qstab:UnknownOption','Unknown option: %s.',name);
    end
end
if ~localTolerance(pairTol) || ~localTolerance(classTol)
    error('qstab:InvalidTolerance','Tolerances must be finite positive real scalars.');
end
[A,B,ci]=qstab.companion(P,'Grade',grade,'Variable',variable);
count=ci.quaternionOrder*ci.grade;
lambda=nan(count,1);
info=struct('status','undetermined','isValid',false,'message','', ...
    'fullSpectrum',zeros(0,1),'pairedIndices',nan(count,2), ...
    'pairingResidual',inf(count,1),'pairingTolerance',pairTol, ...
    'pairingValid',false,'embeddingResidual',ci.embeddingResidual, ...
    'eigenResidual',NaN,'classRepresentatives',zeros(0,1), ...
    'classMultiplicities',zeros(0,1),'classIndex',zeros(0,1), ...
    'classTolerance',classTol,'algebraicCount',count, ...
    'finiteCount',0,'infiniteCount',0,'variable',variable, ...
    'grade',ci.grade,'actualDegree',ci.actualDegree,'companion',ci, ...
    'regularityWitness',[],'spectrumComputed',false, ...
    'assurance','floating-point');
if ci.isZero
    info.status='assumptionFailed';
    info.message='The zero polynomial has no finite algebraic zero list.'; return
end
if ci.grade==0
    lambda=zeros(0,1);
    if ci.constantRcond>100*eps
        info.status='ok'; info.isValid=true; info.pairingValid=true;
        info.eigenResidual=0; info.message='A nonsingular constant polynomial has no zeros.';
        info.regularityWitness=struct('variableValue',0,'reciprocalCondition',ci.constantRcond);
    else
        info.status='assumptionFailed';
        info.message='The constant polynomial matrix is singular or numerically unresolved.';
    end
    return
end
if ci.leadingRcond>100*eps
    info.regularityWitness=struct('kind','leading-coefficient', ...
        'variableValue',Inf,'reciprocalCondition',ci.leadingRcond);
elseif ci.constantRcond>100*eps
    info.regularityWitness=struct('kind','constant-coefficient', ...
        'variableValue',0,'reciprocalCondition',ci.constantRcond);
else
    probe=[0,1,-1,1i,-1i,2,0.5,2i]; best=0; bestShift=NaN;
    for k=1:numel(probe)
        rc=rcond(A-probe(k)*B);
        if rc>best, best=rc; bestShift=probe(k); end
        if rc>100*eps, break; end
    end
    if best<=100*eps
        info.message='No numerically nonsingular pencil value was found; regularity is unresolved.'; return
    end
    info.regularityWitness=struct('kind','pencil-value', ...
        'variableValue',bestShift,'reciprocalCondition',best);
end
[V,D]=builtin('eig',A,B);
full=diag(D);
[lambda,pairs,pairResidual,pairValid,message]=localPair(full,pairTol);
finite=~isnan(full) & ~isinf(full);
if any(finite)
    Vf=V(:,finite); lf=full(finite); VD=bsxfun(@times,Vf,lf.');
    residual=norm(A*Vf-B*VD,'fro')/max(realmin, ...
        norm(A,'fro')*norm(Vf,'fro')+norm(B,'fro')*norm(VD,'fro'));
else
    residual=0;
end
infinite=isinf(full);
if any(infinite)
    Vi=V(:,infinite);
    infResidual=norm(B*Vi,'fro')/max(realmin,norm(B,'fro')*norm(Vi,'fro'));
    residual=max(residual,infResidual);
end
[reps,mult,index]=localGroups(lambda,classTol);
info.fullSpectrum=full; info.pairedIndices=pairs;
info.pairingResidual=pairResidual; info.pairingValid=pairValid;
info.eigenResidual=residual; info.classRepresentatives=reps;
info.classMultiplicities=mult; info.classIndex=index;
info.finiteCount=sum(~isnan(lambda) & ~isinf(lambda));
info.infiniteCount=sum(isinf(lambda)); info.spectrumComputed=true;
info.isValid=pairValid && ~isnan(residual) && ~isinf(residual) && residual<=max(pairTol,1e3*eps);
info.message=message;
if info.isValid
    info.status='ok';
elseif pairValid
    info.message='The generalized eigenvector residual is too large or nonfinite.';
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
