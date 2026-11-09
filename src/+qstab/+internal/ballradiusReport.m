function R = ballradiusReport(P,weights,varargin)
%BALLRADIUSREPORT Internal implementation returning the full report.
ip=inputParser;
addParameter(ip,'Exponent',2);addParameter(ip,'Variable','s',@localText);
addParameter(ip,'Grade',[]);
addParameter(ip,'AbsTol',1e-5,@localPositive);
addParameter(ip,'RelTol',1e-4,@localPositive);
addParameter(ip,'Tolerance',1e-10,@(x)localPositive(x)&&x<1);
addParameter(ip,'MaxCells',8191,@(x)localInteger(x)&&x>=1);
addParameter(ip,'MaxDepth',40,@(x)localInteger(x)&&x>=0);
addParameter(ip,'MaxTime',30,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&~isnan(x)&&x>=0);
addParameter(ip,'StoreCells',true,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});o=ip.Results;variable=lower(char(o.Variable));
if ~any(strcmp(variable,{'s','z','zinv'}))
    error('qstab:ballradius:Variable','Variable must be s, z or zinv.');
end
inputModel=qstab.internal.ballmodel(P,weights,o.Exponent,o.Grade);
M=inputModel;
if strcmp(variable,'zinv')
    M=qstab.internal.ballmodel(localParts(fliplr(inputModel.C)), ...
        fliplr(inputModel.weights),o.Exponent,inputModel.grade);
end
n=M.grade;isCT=strcmp(variable,'s');
if isCT,testVariable='s';domainEnd=1;else,testVariable='z';domainEnd=pi;end
R=struct('schema','qstab.ballradius.1','status','assumptionFailed', ...
    'method','weighted-coefficient-ball-boundary-enclosure', ...
    'variable',variable,'testVariable',testVariable,'grade',n, ...
    'weights',inputModel.weights,'exponent',M.exponent,'validated',false, ...
    'assurance','floating-point bounds, not directed-rounding validation', ...
    'lower',NaN,'upper',NaN,'estimate',NaN,'gap',NaN,'degreeLossRadius',Inf, ...
    'boundaryParameter',NaN,'boundaryArgument',NaN,'limitingEvent','none', ...
    'witness',[],'message','','nominal',struct(),'diagnostics',struct());
timer=tic;
nominal=qstab.internal.isstableReport(localParts(M.C),'Variable',testVariable, ...
    'Grade',n,'Tolerance',o.Tolerance);
R.nominal=nominal;
if ~strcmp(nominal.status,'stable')
    R.message='The nominal polynomial must be numerically verified strictly stable at the declared grade.';
    return
end
if ~M.hasPerturbation
    R.status='unbounded';R.lower=Inf;R.upper=Inf;R.estimate=Inf;R.gap=0;
    R.message='Every coefficient is fixed; the verified nominal member is the whole family.';
    R.diagnostics=struct('processedCells',0,'pointEvaluations',0,'completeDomainCovered',true, ...
        'elapsedSeconds',toc(timer));return
end
% Existence of a zero leading coefficient is a separate event, in testVariable.
upperBound=Inf;bestWitness=[];bestParameter=NaN;bestArgument=NaN;
degreeWitness=[];
if M.weights(end)>0
    degreeRadius=localScale(norm(M.Cscaled(:,end))/M.weightsScaled(end),M.logScaleRatio);
    if isfinite(degreeRadius)&&degreeRadius>0
        R.degreeLossRadius=degreeRadius;
        D=zeros(size(M.C));D(:,end)=-M.C(:,end);Q=M.C;Q(:,end)=0;
        degreeWitness=struct('kind','degreeLoss','valid',true, ...
            'coefficientComponents',Q,'perturbationComponents',D,'polynomial',localParts(Q), ...
            'weightedNorm',degreeRadius,'evaluationResidual',NaN, ...
            'boundaryComponents',[],'classRepresentative',NaN, ...
            'weights',M.weights,'exponent',M.exponent,'grade',n,'validated',false);
        upperBound=degreeRadius;bestWitness=degreeWitness;
        if isCT,bestParameter=1;bestArgument=Inf;end
    end
end
R.degreeLossWitness=localOriginalWitness(degreeWitness,variable,inputModel);
if n==0
    if isempty(degreeWitness)
        R.status='undetermined';R.lower=0;R.upper=Inf;R.gap=Inf;
        R.message='The constant-family radius is outside the representable range.';return
    end
    R.status='converged';R.lower=R.degreeLossRadius;R.upper=R.degreeLossRadius;
    R.estimate=R.upper;R.gap=0;R.witness=R.degreeLossWitness;R.limitingEvent='degreeLoss';
    R.message='Exact-formula constant-family radius; limiting member is the zero polynomial.';
    R.diagnostics=struct('processedCells',0,'pointEvaluations',0,'completeDomainCovered',true, ...
        'elapsedSeconds',toc(timer));return
end
if n>128
    error('qstab:ballradius:DegreeLimit','This bounded solver currently supports declared grades up to 128.');
end
count=0;processed=0;lowerBound=0;reason='Cell or time budget reached.';state='bounded';
if o.MaxTime==0
    R=localFinish(R,0,upperBound,bestWitness,bestParameter,bestArgument, ...
        'bounded','Zero time budget: the full boundary remains unresolved.',variable,inputModel);
    R.diagnostics=struct('processedCells',0,'pointEvaluations',0, ...
        'parameterDomain',[0 domainEnd],'cells',[0 domainEnd 0 0], ...
        'cellColumns',{{'left','right','depth','lower'}}, ...
        'completeDomainCovered',true,'elapsedSeconds',toc(timer));return
end
% Endpoints are explicit, including +1/-1 for DT, w=0 and infinity for CT.
if isCT,seeds=[0 .5];else,seeds=[0 pi/2 pi];end
for t=seeds
    [upperBound,bestWitness,bestParameter,bestArgument]=localUpdatePoint( ...
        M,isCT,t,o.Tolerance,upperBound,bestWitness,bestParameter,bestArgument);
    count=count+1;
end
if isCT
    binom=ones(1,n+1);
    for j=1:n,binom(j+1)=binom(j)*(n-j+1)/j;end
    controls=bsxfun(@times,M.Cscaled,(1i.^(0:n))./binom);
else
    controls=[];
end
shortcut=struct('attempted',false,'accepted',false,'processedCells',0);
if isCT&&isfinite(upperBound)
    target=max(0,upperBound-.8*(o.AbsTol+o.RelTol*upperBound));
    shortcut=qstab.internal.ballpositivity(M,target, ...
        min(1023,floor((o.MaxCells-1)/2)),o.MaxDepth,max(0,o.MaxTime-toc(timer)));
    processed=shortcut.processedCells;
    if shortcut.accepted
        R=localFinish(R,target,upperBound,bestWitness,bestParameter,bestArgument, ...
            'converged','A guarded matrix Bernstein inequality supplies the global lower radius.',variable,inputModel);
        R.diagnostics=struct('processedCells',processed,'leafCount',size(shortcut.acceptedCells,1), ...
            'pointEvaluations',count,'parameterDomain',[0 1], ...
            'parameterMeaning','t=omega/(1+omega); t=1 is the infinite-frequency endpoint', ...
            'lowerBoundMethod','Bernstein matrix positivity at a candidate radius', ...
            'completeDomainCovered',true,'positivity',shortcut, ...
            'elapsedSeconds',toc(timer),'validated',false);
        return
    end
end
% Cell evaluation is a file-local function with explicit inputs, NOT a
% nested function sharing the subdivision variables of this search.
[upperBound,bestWitness,bestParameter,bestArgument]=localUpdatePoint( ...
    M,isCT,domainEnd/2,o.Tolerance,upperBound,bestWitness,bestParameter,bestArgument);
count=count+1;processed=processed+1;
cells=localCell(M,isCT,0,domainEnd,0,controls);
while true
    lowers=[cells.lower];lowerBound=min(lowers);
    % The separate degree-loss member is also a possible limiting event.
    lowerBound=min(lowerBound,R.degreeLossRadius);
    if lowerBound>upperBound
        if lowerBound-upperBound>2048*eps*max(upperBound,realmin)
            state='undetermined';reason='Inconsistent numerical lower and upper bounds.';break
        end
        lowerBound=max(0,upperBound-2048*eps*max(upperBound,realmin));
    end
    if isfinite(upperBound)&&upperBound-lowerBound<=o.AbsTol+o.RelTol*upperBound
        state='converged';reason='The global numerical radius bracket meets the requested tolerance.';break
    end
    if processed+2>o.MaxCells||toc(timer)>=o.MaxTime,break;end
    [~,idx]=min(lowers);old=cells(idx);splitPoint=old.left+(old.right-old.left)/2;
    if old.depth>=o.MaxDepth
        reason='The lowest-bound cell reached MaxDepth; it remains an unresolved obligation.';break
    end
    if splitPoint==old.left||splitPoint==old.right
        reason='Subdivision reached floating-point resolution.';break
    end
    if isCT,[cl,cr]=splitControls(old.controls);else,cl=[];cr=[];end
    leftMidpoint=old.left+(splitPoint-old.left)/2;
    [upperBound,bestWitness,bestParameter,bestArgument]=localUpdatePoint( ...
        M,isCT,leftMidpoint,o.Tolerance,upperBound,bestWitness,bestParameter,bestArgument);
    count=count+1;processed=processed+1;
    left=localCell(M,isCT,old.left,splitPoint,old.depth+1,cl);
    rightMidpoint=splitPoint+(old.right-splitPoint)/2;
    [upperBound,bestWitness,bestParameter,bestArgument]=localUpdatePoint( ...
        M,isCT,rightMidpoint,o.Tolerance,upperBound,bestWitness,bestParameter,bestArgument);
    count=count+1;processed=processed+1;
    right=localCell(M,isCT,splitPoint,old.right,old.depth+1,cr);
    cells(idx)=left;cells(end+1)=right; %#ok<AGROW>
end
if ~isfinite(upperBound),state='undetermined';reason='No representable finite witness upper bound was reconstructed.';end
% Never attach a full-domain claim to a malformed partition.  This guard is
% independent of StoreCells and does not relax any numerical tolerance.
partitionVerified=localPartition(cells,domainEnd);
if ~partitionVerified
    lowerBound=0;state='undetermined';
    reason='The final cells do not form an exact partition of the full boundary.';
end
R=localFinish(R,lowerBound,upperBound,bestWitness,bestParameter,bestArgument, ...
    state,reason,variable,inputModel);
if o.StoreCells
    table=[[cells.left].',[cells.right].',[cells.depth].',[cells.lower].', ...
        [cells.minimumModulus].',[cells.variation].',[cells.envelopeUpper].'];
else
    table=zeros(0,7);
end
R.diagnostics=struct('processedCells',processed,'leafCount',numel(cells), ...
    'pointEvaluations',count,'parameterDomain',[0 domainEnd], ...
    'parameterMeaning','','cells',table, ...
    'cellColumns',{{'left','right','depth','lower','scaledSigmaMin','scaledVariation','scaledEnvelopeUpper'}}, ...
    'cellsStored',o.StoreCells,'completeDomainCovered',partitionVerified, ...
    'partitionVerified',partitionVerified,'positivity',shortcut, ...
    'AbsTol',o.AbsTol,'RelTol',o.RelTol,'MaxCells',o.MaxCells,'MaxDepth',o.MaxDepth, ...
    'MaxTime',o.MaxTime,'elapsedSeconds',toc(timer),'validated',false);
if isCT
    R.diagnostics.parameterMeaning='t=omega/(1+omega); t=1 is the infinite-frequency endpoint';
    R.diagnostics.lowerBoundMethod='homogeneous Bernstein matrix variation / Holder envelope upper bound';
else
    R.diagnostics.parameterMeaning='theta in [0,pi], tested forward z=exp(1i*theta)';
    R.diagnostics.lowerBoundMethod='unit-circle matrix Taylor enclosure with local unit-modulus phase';
end

end
function [upperBound,bestWitness,bestParameter,bestArgument]=localUpdatePoint( ...
    M,isCT,t,tolerance,upperBound,bestWitness,bestParameter,bestArgument)
% Explicit state in/out prevents accidental capture of optimizer variables.
if isCT,z=1i*t/(1-t);
elseif t==0,z=1;
elseif t==pi,z=-1;
else,z=exp(1i*t);end
E=qstab.internal.ballpoint(M,z,true,tolerance);
if E.valid&&~isempty(E.witness)&&E.witness.valid&&isfinite(E.upper)&&E.upper<upperBound
    upperBound=E.upper;bestWitness=E.witness;bestParameter=t;bestArgument=z;
end
end
function cellOut=localCell(M,isCT,l,r,depth,Cctrl)
% Pure local cell calculation: no shared workspace or side effects.
n=M.grade;mid=l+(r-l)/2;half=(r-l)/2;
    if isCT
        v=midControls(Cctrl);B=localMatrix(v);s=svd(B);m=s(end);
        variation=0;amp=0;
        for ii=1:n+1
            variation=max(variation,norm(localMatrix(Cctrl(:,ii)-v),2));
            amp=max(amp,norm(localMatrix(Cctrl(:,ii)),2));
        end
        maxPowers=zeros(1,n+1);
        for kk=0:n
            peak=min(r,max(l,kk/n));
            maxPowers(kk+1)=peak^kk*(1-peak)^(n-kk);
        end
        den=qstab.internal.ballnorm(M.weightsScaled.*maxPowers,M.dualExponent);
        err=256*eps*(n+1)*(depth+1)*max(amp,realmin);
    else
        kvec=0:n;phase=exp(1i*kvec*mid);
        v=M.Cscaled*phase.';dv=M.Cscaled*(1i*kvec.*phase).';
        B=localMatrix(v);D=localMatrix(dv);s=svd(B);m=s(end);
        norms=zeros(1,n+1);
        for ii=1:n+1,norms(ii)=norm(M.Cscaled(:,ii));end
        phases=[0 n/2];denB=norm(B,'fro')^2;
        if denB>0
            alpha=imag(sum(conj(B(:)).*D(:)))/denB;
            if isfinite(alpha),phases(end+1)=alpha;end
        end
        variation=Inf;
        for alpha=phases
            first=norm(D-1i*alpha*B,2);
            second=sum((kvec-alpha).^2.*norms);
            variation=min(variation,half*first+half^2*second/2);
        end
        amp=sum(norms);den=qstab.internal.ballnorm(M.weightsScaled,M.dualExponent);
        err=256*eps*(n+1)*(depth+1)*max(amp+variation,realmin);
    end
    den=den*(1+256*eps*(n+1)*(depth+1));
    numerator=max(0,m-variation-err);lb=0;
    if isfinite(numerator)&&isfinite(den)&&den>0
        lb=localScale(numerator/den,M.logScaleRatio);
        if ~isfinite(lb),lb=0;end
    end
    cellOut=struct('left',l,'right',r,'depth',depth,'lower',lb, ...
        'minimumModulus',m,'variation',variation,'envelopeUpper',den,'controls',Cctrl);
end
function ok=localPartition(cells,domainEnd)
% Shared endpoints are copied, so exact equality is intentional here.
T=sortrows([[cells.left].',[cells.right].'],1);
ok=~isempty(T)&&all(isfinite(T(:)))&&T(1,1)==0&&T(end,2)==domainEnd ...
    &&all(T(:,2)>T(:,1))&&all(T(1:end-1,2)==T(2:end,1));
end
function R=localFinish(R,lb,ub,W,t,z,state,reason,variable,inputModel)
R.status=state;R.lower=lb;R.upper=ub;R.estimate=ub;R.gap=ub-lb;
R.witness=localOriginalWitness(W,variable,inputModel);
R.boundaryParameter=t;R.boundaryArgument=z;R.message=reason;
if ~isempty(W),R.limitingEvent=W.kind;end
if strcmp(variable,'zinv')&&~isempty(W)&&strcmp(W.kind,'finiteZero')
    R.boundaryArgument=1/z;
end
end
function W=localOriginalWitness(W,variable,inputModel)
if isempty(W),return;end
W.variable=variable;W.testedCoefficientComponents=W.coefficientComponents;
W.testedPerturbationComponents=W.perturbationComponents;
W.testedBoundaryComponents=W.boundaryComponents;
W.testedClassRepresentative=W.classRepresentative;
if strcmp(variable,'zinv')
    W.coefficientComponents=fliplr(W.coefficientComponents);
    W.perturbationComponents=fliplr(W.perturbationComponents);
    W.polynomial=localParts(W.coefficientComponents);W.weights=inputModel.weights;
    if strcmp(W.kind,'finiteZero')
        q=W.boundaryComponents;
        W.boundaryComponents=[q(1);-q(2:4)]/sum(q.^2);
        W.classRepresentative=complex(W.boundaryComponents(1),norm(W.boundaryComponents(2:4)));
    end
end
if strcmp(W.kind,'degreeLoss')
    if strcmp(variable,'zinv'),W.eventDescription='Loss of the current-sample coefficient (forward leading coefficient), not a finite delay zero.';
    else,W.eventDescription='Loss of the required leading coefficient, not a finite boundary-zero claim.';end
else
    W.eventDescription='Numerically reconstructed algebraic boundary-zero member of the original input family.';
end
end
function [L,R]=splitControls(C)
n=size(C,2);L=zeros(size(C));R=L;D=C;L(:,1)=D(:,1);R(:,n)=D(:,n);
for k=2:n
    D=(D(:,1:end-1)+D(:,2:end))/2;L(:,k)=D(:,1);R(:,n-k+1)=D(:,end);
end
end
function v=midControls(C)
while size(C,2)>1,C=(C(:,1:end-1)+C(:,2:end))/2;end
v=C(:,1);
end
function B=localMatrix(v)
B=[v(1)+1i*v(2),v(3)+1i*v(4);-v(3)+1i*v(4),v(1)-1i*v(2)];
end
function y=localScale(x,logScale)
if x==0,y=0;elseif isinf(x),y=Inf;else,y=exp(log(x)+logScale);end
end
function P=localParts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
function tf=localText(x)
tf=(ischar(x)&&isrow(x))||(isstring(x)&&isscalar(x));
end
function tf=localPositive(x)
tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0;
end
function tf=localInteger(x)
tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x==fix(x);
end
