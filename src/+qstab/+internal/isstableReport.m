function result = isstableReport(X,varargin)
%ISSTABLEREPORT Internal implementation returning the full report.

o=localOptions(varargin{:});
[result,boundaryInput]=localCompute(X,o);
% Keep the original, unscaled coefficients and all original diagnostics.
% A tolerance-band zero or a zero computed after normalisation is not an
% exact certificate. The structural fallback needs no eigensolver.
if strcmp(result.status,'undetermined') && ~isempty(boundaryInput) && ...
        ~isfield(result.diagnostics,'errorIdentifier')
    W=qstab.internal.exactBoundary(o.type,boundaryInput,o.domain);
    if ~isempty(W)
        result=localDecision(result,'notStable',W.reason);
        result.diagnostics.boundaryWitness=W;
    end
end
end

function [result,boundaryInput]=localCompute(X,o)
boundaryInput=[];
result=localResult(o);
try
    switch o.type
        case 'matrix'
            [M,mi]=qstab.chi(X);
            if isempty(M) || size(M,1)~=size(M,2)
                error('qstab:NotSquareMatrix','A nonempty square matrix is required.');
            end
            if ~isempty(o.grade)
                error('qstab:InvalidGrade','Grade applies to polynomial inputs only.');
            end
            result.diagnostics.embedding=mi; boundaryInput=X;
            if strcmp(o.method,'spectral')
                [lambda,si]=qstab.eig(X,'Tolerance',max(o.tolerance,1e3*eps));
                result.spectrumComputed=true; result.method='spectral';
                result.diagnostics.spectrum=si; result.diagnostics.representatives=lambda;
                if ~si.isValid
                    result.message=si.message; return
                end
                [status,msg,gap]=localSpectrum(lambda,o.domain,o.tolerance);
                result=localDecision(result,status,msg); result.diagnostics.margin=gap;
            elseif ismember(o.method,{'auto','lyapunov'})
                [status,msg,di]=localMatrix(M,o);
                result=localDecision(result,status,msg); result.method='lyapunov';
                result.diagnostics.certificate=di;
            else
                error('qstab:InvalidMethod','Matrix inputs require auto, lyapunov, or spectral.');
            end
        case 'polynomial'
            [C,pi]=qstab.polycomponents(X,'Grade',o.grade);
            result.grade=pi.grade; result.actualDegree=pi.actualDegree;
            result.diagnostics.input=pi;
            if pi.isZero
                error('qstab:ZeroPolynomial','The zero polynomial has no strict stability test.');
            end
            if strcmp(o.variable,'zinv')
                C=fliplr(C); result.diagnostics.reversedAtGrade=pi.grade;
                result.diagnostics.testVariable='z';
            else
                result.diagnostics.testVariable=o.variable;
            end
            if all(C(:,end)==0)
                error('qstab:LeadingCoefficientMissing', ...
                    'The leading coefficient in the forward variable vanishes at the declared grade.');
            end
            if (strcmp(o.domain,'continuous') && strcmp(o.method,'schur')) || ...
               (strcmp(o.domain,'discrete') && strcmp(o.method,'routh'))
                error('qstab:InvalidMethod','Method and stability domain do not agree.');
            end
            boundaryInput=C; % original forward coefficients, before scaling
            scale=max(abs(C(:))); C=C/scale;
            result.diagnostics.inputScale=scale;
            if norm(C(:,end))<=o.tolerance
                result.message='The leading coefficient is too small relative to the supplied coefficients.'; return
            end
            P=localComponentStruct(C);
            if strcmp(o.method,'spectral')
                [lambda,si]=qstab.roots(P,'Grade',pi.grade, ...
                    'Variable',result.diagnostics.testVariable, ...
                    'Tolerance',max(o.tolerance,1e3*eps));
                result.spectrumComputed=si.spectrumComputed; result.method='spectral';
                result.diagnostics.spectrum=si; result.diagnostics.representatives=lambda;
                if ~si.isValid
                    result.message=si.message; return
                end
                [status,msg,gap]=localSpectrum(lambda,o.domain,o.tolerance);
                result=localDecision(result,status,msg); result.diagnostics.margin=gap;
            elseif strcmp(o.method,'lyapunov')
                M=qstab.companion(P,'Grade',pi.grade);
                if isempty(M)
                    result=localDecision(result,'stable','A nonzero constant polynomial has no unstable modes.');
                    result.method='constant'; return
                end
                [status,msg,di]=localMatrix(M,o);
                result=localDecision(result,status,msg); result.method='companion-lyapunov';
                result.diagnostics.certificate=di;
            else
                [np,ni]=qstab.normpoly(P,'Grade',pi.grade);
                result.diagnostics.normPolynomial=np;
                result.diagnostics.normPolynomialInfo=ni;
                if strcmp(o.domain,'continuous')
                    if strcmp(o.method,'schur')
                        error('qstab:InvalidMethod','Schur recursion requires z or zinv.');
                    end
                    [status,msg,di]=localRouth(np,o.tolerance); result.method='routh';
                else
                    if strcmp(o.method,'routh')
                        error('qstab:InvalidMethod','Routh recursion requires s.');
                    end
                    [status,msg,di]=localSchur(np,o.tolerance); result.method='schur';
                end
                result=localDecision(result,status,msg); result.diagnostics.test=di;
            end
        case 'polynomialmatrix'
            if ~ismember(o.method,{'auto','lyapunov','spectral'})
                error('qstab:InvalidMethod','Polynomial matrices require auto, lyapunov, or spectral.');
            end
            P=X;
            if strcmp(o.variable,'zinv')
                P=localReversePolynomialMatrix(X,o.grade);
                result.diagnostics.testVariable='z';
            else
                result.diagnostics.testVariable=o.variable;
            end
            [Ap,Bp,ci]=qstab.companion(P,'Grade',o.grade);
            result.grade=ci.grade; result.actualDegree=ci.actualDegree;
            result.diagnostics.input=ci;
            if ci.leadingRcond<=100*eps
                error('qstab:SingularLeadingCoefficient', ...
                    'This nominal test requires a nonsingular forward leading coefficient.');
            end
            if ci.grade==0
                result=localDecision(result,'stable','A nonsingular constant polynomial matrix has no zeros.');
                result.method='constant'; return
            end
            % The top companion row and the leading B block are signed
            % copies of the original complex adjoints, unlike Bp\Ap.
            boundaryInput=struct('A',Ap,'B',Bp,'n',ci.quaternionOrder,'grade',ci.grade);
            if strcmp(o.method,'spectral')
                [lambda,si]=qstab.roots(P,'Grade',ci.grade, ...
                    'Variable',result.diagnostics.testVariable, ...
                    'Tolerance',max(o.tolerance,1e3*eps));
                result.spectrumComputed=si.spectrumComputed; result.method='spectral';
                result.diagnostics.spectrum=si; result.diagnostics.representatives=lambda;
                if ~si.isValid, result.message=si.message; return; end
                [status,msg,gap]=localSpectrum(lambda,o.domain,o.tolerance);
                result=localDecision(result,status,msg); result.diagnostics.margin=gap;
            else
                M=Bp\Ap;
                [status,msg,di]=localMatrix(M,o);
                result=localDecision(result,status,msg); result.method='companion-lyapunov';
                result.diagnostics.certificate=di;
            end
    end
catch ME
    if strncmp(ME.identifier,'qstab:',6)
        result=localDecision(result,'assumptionFailed',ME.message);
    else
        result=localDecision(result,'undetermined',ME.message);
    end
    result.diagnostics.errorIdentifier=ME.identifier;
end
end

function o=localOptions(varargin)
o=struct('type','polynomial','variable','s','domain','', ...
    'grade',[],'method','auto','tolerance',1e-10,'maxKroneckerOrder',24);
hasVariable=false;
if mod(numel(varargin),2)~=0
    error('qstab:InvalidOptions','Options must be name/value pairs.');
end
for k=1:2:numel(varargin)
    name=lower(char(varargin{k})); value=varargin{k+1};
    switch name
        case 'type', o.type=lower(char(value));
        case 'variable', o.variable=lower(char(value)); hasVariable=true;
        case 'domain', o.domain=lower(char(value));
        case 'grade', o.grade=value;
        case 'method', o.method=lower(char(value));
        case 'tolerance', o.tolerance=value;
        case 'maxkroneckerorder', o.maxKroneckerOrder=value;
        otherwise, error('qstab:UnknownOption','Unknown option: %s.',name);
    end
end
if ~ismember(o.type,{'matrix','polynomial','polynomialmatrix'})
    error('qstab:InvalidType','Type must be matrix, polynomial, or polynomialMatrix.');
end
if ~ismember(o.variable,{'s','z','zinv'})
    error('qstab:InvalidVariable','Variable must be s, z, or zinv.');
end
if isempty(o.domain)
    if strcmp(o.variable,'s'), o.domain='continuous'; else, o.domain='discrete'; end
elseif ~ismember(o.domain,{'continuous','discrete'})
    error('qstab:InvalidDomain','Domain must be continuous or discrete.');
elseif ~hasVariable
    if strcmp(o.domain,'continuous'), o.variable='s'; else, o.variable='z'; end
elseif strcmp(o.domain,'continuous')~=strcmp(o.variable,'s')
    error('qstab:InconsistentDomain','Domain and Variable describe different stability regions.');
end
if strcmp(o.type,'matrix') && strcmp(o.variable,'zinv')
    error('qstab:InvalidMatrixVariable', ...
        'A state transition matrix is tested in forward time; use Domain discrete.');
end
if ismember(o.method,{'eig','roots','eigenvalues'}), o.method='spectral'; end
if ~ismember(o.method,{'auto','coefficient','routh','schur','lyapunov','spectral'})
    error('qstab:InvalidMethod','Unknown stability method.');
end
if ~(isnumeric(o.tolerance) && isreal(o.tolerance) && isscalar(o.tolerance) && ...
        ~isnan(o.tolerance) && ~isinf(o.tolerance) && o.tolerance>0 && o.tolerance<1)
    error('qstab:InvalidTolerance','Tolerance must lie strictly between zero and one.');
end
if ~(isnumeric(o.maxKroneckerOrder) && isreal(o.maxKroneckerOrder) && ...
        isscalar(o.maxKroneckerOrder) && ~isnan(o.maxKroneckerOrder) && ...
        ~isinf(o.maxKroneckerOrder) && o.maxKroneckerOrder>=1 && ...
        o.maxKroneckerOrder==fix(o.maxKroneckerOrder))
    error('qstab:InvalidOptions','MaxKroneckerOrder must be a positive integer.');
end
end

function r=localResult(o)
r=struct('status','undetermined','isStable',[],'message','', ...
    'type',o.type,'variable',o.variable,'domain',o.domain, ...
    'method',o.method,'grade',[],'actualDegree',[], ...
    'spectrumComputed',false,'assurance','floating-point', ...
    'tolerance',o.tolerance,'diagnostics',struct());
end

function r=localDecision(r,status,message)
r.status=status; r.message=message;
if strcmp(status,'stable'), r.isStable=true;
elseif strcmp(status,'notStable'), r.isStable=false;
else, r.isStable=[];
end
end

function P=localComponentStruct(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end

function [status,msg,di]=localRouth(p,tol)
% Keep every resolved table decision. Only an unresolved quartic receives
% the one-sided algebraic fallback below; it can never assert stability.
[status,msg,di]=localRouthTable(p,tol);
if strcmp(status,'undetermined') && numel(p)==5
    [unstable,witness]=localQuarticHurwitzWitness(p,tol);
    di.hurwitzFallback=witness;
    if unstable
        status='notStable';
        msg=['A separated negative third Hurwitz determinant proves a ', ...
            'right-half-plane root despite the unresolved Routh table.'];
    end
end
end

function [unstable,W]=localQuarticHurwitzWitness(p,tol)
% Descending coefficients [a4 a3 a2 a1 a0] after positive scaling:
% Delta3 = a3*a2*a1 - a4*a1^2 - a3^2*a0.
% For roots in the CLOSED left half-plane every Hurwitz determinant is
% nonnegative: shift p(s) to p(s+epsilon), epsilon>0, then take the limit.
% A safely negative Delta3 therefore proves a strictly right-half-plane
% root. Delta3 near zero proves nothing and is deliberately not replaced
% by a sign, a small pivot or an eigenvalue-based answer.
p=p(:).';scale=max(abs(p));
W=struct('name','quarticThirdHurwitzDeterminant','order',3, ...
    'coefficientScale',scale,'coefficientsDescending',[], ...
    'signedTerms',[],'normalizedDeterminant',NaN, ...
    'errorAllowance',NaN,'negativeGap',NaN, ...
    'provesRightHalfPlaneRoot',false,'validated',false);
unstable=false;
if scale==0 || any(~isfinite(p)),return;end
c=fliplr(p/scale);if c(1)<0,c=-c;end
terms=[c(2)*c(3)*c(4),-c(1)*c(4)*c(4),-c(2)*c(2)*c(5)];
delta=(terms(1)+terms(2))+terms(3);
% Ordinary double precision: a conservative sign guard, NOT an
% outward-rounded interval proof. Normalization makes |c_j|<=1.
allowance=(tol+64*eps)*max(1,sum(abs(terms)));
W.coefficientsDescending=c;W.signedTerms=terms;
W.normalizedDeterminant=delta;W.errorAllowance=allowance;
W.negativeGap=-delta-allowance;
unstable=all(isfinite([delta,allowance])) && delta < -allowance;
W.provesRightHalfPlaneRoot=unstable;
end

function [status,msg,di]=localRouthTable(p,tol)
p=p(:).'; p=p/max(abs(p)); d=numel(p)-1;
di=struct('degree',d,'firstColumn',zeros(0,1),'criticalRow',[]);
status='undetermined'; msg='The Routh test is numerically unresolved.';
if d==0, status='stable'; msg='A nonzero constant polynomial has no roots.'; return; end
p=fliplr(p);
if p(1)<0, p=-p; end
if any(p < -tol)
    status='notStable'; msg='A negative normalized coefficient excludes closed-left-half-plane stability.'; return
end
if any(abs(p)<=tol)
    msg='A normalized coefficient is zero or near zero; no boundary decision is made.'; return
end
cols=ceil((d+1)/2); R=zeros(d+1,cols);
R(1,1:numel(p(1:2:end)))=p(1:2:end);
R(2,1:numel(p(2:2:end)))=p(2:2:end);
R(1,:)=R(1,:)/max(abs(R(1,:)));
R(2,:)=R(2,:)/max(abs(R(2,:)));
for k=3:d+1
    pivot=R(k-1,1);
    if abs(pivot)<=tol
        di.firstColumn=R(1:k-1,1); di.criticalRow=k-1;
        msg='A Routh pivot is zero or near zero; no epsilon substitution is used.'; return
    end
    for j=1:cols-1
        R(k,j)=(pivot*R(k-2,j+1)-R(k-2,1)*R(k-1,j+1))/pivot;
    end
    rowScale=max(abs(R(k,:)));
    if rowScale<=tol || isnan(rowScale) || isinf(rowScale)
        di.firstColumn=R(1:k,1); di.criticalRow=k;
        msg='A Routh row vanishes or is numerically unresolved.'; return
    end
    R(k,:)=R(k,:)/rowScale;
    if R(k,1)<-tol
        status='notStable'; msg='A resolved Routh sign change detects a right-half-plane root.';
        di.firstColumn=R(1:k,1); di.criticalRow=k; return
    end
end
di.firstColumn=R(:,1);
if all(R(:,1)>tol)
    status='stable'; msg='All normalized Routh pivots are positive and separated from zero.';
else
    msg='A final Routh pivot lies in the unresolved tolerance band.';
end
end

function [status,msg,di]=localSchur(p,tol)
p=p(:).'; degree=numel(p)-1; gaps=zeros(degree,1);
status='undetermined'; msg='The Schur recursion is numerically unresolved.';
di=struct('degree',degree,'reflectionGaps',gaps,'completedSteps',0);
for k=1:degree
    scale=max(abs(p));
    if scale==0 || isnan(scale) || isinf(scale)
        msg='The reduced Schur polynomial vanishes or is nonfinite.'; return
    end
    p=p/scale;
    gap=abs(p(end))-abs(p(1)); gaps(k)=gap;
    di.reflectionGaps=gaps; di.completedSteps=k;
    if gap < -tol
        status='notStable'; msg='A resolved Schur inequality excludes the closed unit disk.'; return
    elseif gap<=tol
        msg='A Schur inequality is on or near equality; the unit-circle case is unresolved.'; return
    end
    p=p(end)*p(2:end)-p(1)*fliplr(p(1:end-1));
end
status='stable'; msg='Every Schur recursion inequality is strictly satisfied with numerical separation.';
end

function [status,msg,di]=localSpectrum(lambda,domain,tol)
di=NaN; status='undetermined'; msg='The computed spectrum is unresolved.';
if any(isnan(lambda) | isinf(lambda)), return; end
if isempty(lambda), status='stable'; msg='There are no finite dynamic modes.'; di=Inf; return; end
if strcmp(domain,'continuous')
    edge=max(real(lambda)); band=tol*max(1,max(abs(lambda))); di=-edge;
else
    edge=max(abs(lambda))-1; band=tol*max(1,max(abs(lambda))); di=-edge;
end
if edge < -band
    status='stable'; msg='All computed representatives are inside the stability region with numerical separation.';
elseif edge>band
    status='notStable'; msg='A computed representative is outside the stability region with numerical separation.';
else
    msg='A computed representative lies on or near the stability boundary.';
end
end

function [status,msg,di]=localMatrix(M,o)
n=size(M,1); I=eye(n); tol=o.tolerance;
di=struct('solver','','P',[],'Q',[],'equationResidual',NaN, ...
    'reciprocalCondition',NaN,'marginLowerBound',NaN,'negativeDirection',[]);
status='undetermined'; msg='A Lyapunov certificate is numerically unresolved.';
if strcmp(o.domain,'continuous')
    Q=-(M'+M);
else
    Q=I-M'*M;
end
[positive,qLower]=localPositive(Q);
if positive
    gap=localGap(qLower,1,o.domain,true);
    if gap>tol*localBoundaryScale(M,o.domain)
        di.solver='identity-certificate'; di.P=I; di.Q=Q;
        di.equationResidual=0; di.marginLowerBound=gap;
        status='stable'; msg='The identity matrix provides a separated Lyapunov/Stein certificate.'; return
    end
end
P=[];
if n<=o.maxKroneckerOrder
    if strcmp(o.domain,'continuous')
        K=kron(I,M')+kron(M.',I);
    else
        K=kron(M.',M')-eye(n*n);
    end
    rc=rcond(K); di.reciprocalCondition=rc; di.solver='kronecker';
    if rc<=100*eps
        msg='The Lyapunov/Stein operator is singular or too ill-conditioned.'; return
    end
    rhs=-I; P=reshape(K\rhs(:),n,n);
else
    if strcmp(o.domain,'continuous'), solver='lyap'; else, solver='dlyap'; end
    if exist(solver,'file')~=0
        try
            P=feval(solver,M',I); di.solver=solver;
        catch ME
            msg=['The optional Lyapunov solver failed: ',ME.message]; return
        end
    else
        msg='The base-MATLAB Kronecker size limit is exceeded and no optional Lyapunov solver is available.'; return
    end
end
P=(P+P')/2;
if any(isnan(P(:)) | isinf(P(:)))
    msg='The computed Lyapunov matrix is nonfinite.'; return
end
if strcmp(o.domain,'continuous'), Q=-(M'*P+P*M); else, Q=P-M'*P*M; end
Q=(Q+Q')/2;
di.P=P; di.Q=Q; di.equationResidual=norm(Q-I,'fro')/max(1,sqrt(n));
[qPositive,qLower]=localPositive(Q);
if ~qPositive
    msg='The recomputed Lyapunov/Stein right-hand side is not reliably positive definite.'; return
end
[pPositive,~]=localPositive(P);
pUpper=norm(P,1);
if pPositive
    gap=localGap(qLower,pUpper,o.domain,true); di.marginLowerBound=gap;
    if gap>tol*localBoundaryScale(M,o.domain)
        status='stable'; msg='A positive Lyapunov/Stein certificate has a separated stability margin.';
    else
        msg='The certificate gives insufficient separation from the stability boundary.';
    end
else
    [negative,v]=localNegativeDirection(P,tol);
    if negative
        gap=localGap(qLower,pUpper,o.domain,false); di.marginLowerBound=gap;
        di.negativeDirection=v;
        if gap>tol*localBoundaryScale(M,o.domain)
            status='notStable'; msg='A resolved negative direction of the Lyapunov matrix proves a separated unstable mode.';
        else
            msg='The Lyapunov inertia indicates loss of stability, but the boundary separation is unresolved.';
        end
    else
        msg='Positive definiteness of the Lyapunov matrix is unresolved.';
    end
end
end

function scale=localBoundaryScale(M,domain)
if strcmp(domain,'continuous'), scale=max(1,norm(M,1)); else, scale=1; end
end

function gap=localGap(qLower,pUpper,domain,isStable)
if pUpper<=0 || isnan(pUpper) || isinf(pUpper), gap=0; return; end
ratio=qLower/pUpper;
if strcmp(domain,'continuous')
    gap=ratio/2;
elseif isStable
    ratio=min(1,max(0,ratio)); gap=ratio/(1+sqrt(1-ratio));
else
    gap=ratio/(sqrt(1+ratio)+1);
end
end

function [tf,lower]=localPositive(H)
H=(H+H')/2; tf=false; lower=0;
if any(isnan(H(:)) | isinf(H(:))) || isempty(H), return; end
scale=norm(H,1);
if scale<=0 || isnan(scale) || isinf(scale), return; end
Hs=H/scale;
[~,flag]=chol(Hs);
if flag~=0 || rcond(Hs)<=100*eps, return; end
invNorm=norm(Hs\eye(size(Hs)),1);
lower=scale/invNorm;
tf=lower>0 && ~isnan(lower) && ~isinf(lower);
end

function [tf,v]=localNegativeDirection(P,tol)
tf=false; v=[]; n=size(P,1); scale=norm(P,1);
if scale<=0, return; end
for k=1:n
    if real(P(k,k)) < -tol*scale
        v=zeros(n,1); v(k)=1; tf=true; return
    end
end
try
    [L,D,perm]=ldl(P,'vector');
catch
    return
end
k=1;
while k<=n
    y=zeros(n,1);
    if k<n && D(k+1,k)~=0
        a=real(D(k,k)); d=real(D(k+1,k+1)); b=D(k,k+1);
        mu=(a+d-hypot(a-d,2*abs(b)))/2;
        if mu<0
            w=[b;mu-a];
            if norm(w)>0, y(k:k+1)=w/norm(w); end
        end
        k=k+2;
    else
        if real(D(k,k))<0, y(k)=1; end
        k=k+1;
    end
    if any(y~=0)
        w=L'\y; candidate=zeros(n,1); candidate(perm)=w;
        candidate=candidate/norm(candidate);
        if real(candidate'*P*candidate)<-tol*scale
            v=candidate; tf=true; return
        end
    end
end
end

function P=localReversePolynomialMatrix(X,grade)
if iscell(X)
    if isempty(X) || ~isvector(X), error('qstab:InvalidPolynomial','Use a coefficient cell vector.'); end
    P=X(:).'; inputGrade=numel(P)-1;
    if isempty(grade), grade=inputGrade; end
    if ~localGradeValid(grade,inputGrade)
        error('qstab:InvalidGrade','Grade must preserve all supplied coefficients.');
    end
    [a,~,~,~]=qstab.components(P{1});
    for k=inputGrade+2:grade+1, P{k}=zeros(size(a)); end
    P=fliplr(P);
else
    [a,b,c,d]=qstab.components(X);
    if isvector(a) || isscalar(a)
        a=reshape(a,1,1,[]); b=reshape(b,1,1,[]);
        c=reshape(c,1,1,[]); d=reshape(d,1,1,[]);
    end
    inputGrade=size(a,3)-1;
    if isempty(grade), grade=inputGrade; end
    if ~localGradeValid(grade,inputGrade)
        error('qstab:InvalidGrade','Grade must preserve all supplied coefficients.');
    end
    if grade>inputGrade
        a(:,:,end+1:grade+1)=0; b(:,:,end+1:grade+1)=0;
        c(:,:,end+1:grade+1)=0; d(:,:,end+1:grade+1)=0;
    end
    P=struct('a',a(:,:,end:-1:1),'b',b(:,:,end:-1:1), ...
        'c',c(:,:,end:-1:1),'d',d(:,:,end:-1:1));
end
end

function tf=localGradeValid(g,inputGrade)
tf=isnumeric(g) && isreal(g) && isscalar(g) && ~isnan(g) && ~isinf(g) && ...
    g>=inputGrade && g==fix(g);
end
