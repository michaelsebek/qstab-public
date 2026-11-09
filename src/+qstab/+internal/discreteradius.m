function R = discreteradius(P0,P1,o)
%DISCRETERADIUS Global ONE quaternion parameter radius for z and zinv.
% This is not the coefficient-ball problem. All coefficients share one t.
side=lower(char(o.Side));variable=lower(char(o.Variable));
if ~any(strcmp(side,{'left','right'})),error('qstab:stabilityradius:Side','Side must be left or right.');end
positive={'AbsTol','RelTol','Tolerance'};
for j=1:numel(positive)
    x=o.(positive{j});
    if ~(isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0)
        error('qstab:stabilityradius:Options','Tolerances must be finite positive scalars.');
    end
end
if o.Tolerance>=1,error('qstab:stabilityradius:Options','Tolerance must be less than one.');end
if ~integer(o.MaxCells,1)||~integer(o.MaxDepth,0)|| ...
        ~(isnumeric(o.MaxTime)&&isreal(o.MaxTime)&&isscalar(o.MaxTime)&&~isnan(o.MaxTime)&&o.MaxTime>=0)|| ...
        ~(islogical(o.StoreCells)&&isscalar(o.StoreCells))
    error('qstab:stabilityradius:Budget','Invalid discrete radius budget or StoreCells.');
end
M=qstab.internal.discretepair(P0,P1,variable,o.Grade);n=M.grade;
if n>64,error('qstab:stabilityradius:DegreeLimit','The discrete radius supports declared grades up to 64.');end
R=struct('schema','qstab.discrete-radius.1','status','assumptionFailed', ...
    'method','one-quaternion-unit-boundary-enclosure','assurance','floating-point-bounds', ...
    'validated',false,'side',side,'variable',variable,'testVariable','z','grade',n, ...
    'lower',NaN,'upper',NaN,'estimate',NaN,'gap',NaN,'frequency',NaN, ...
    'boundaryParameter',NaN,'boundaryArgument',NaN,'limitingEvent','none', ...
    'criticalParameterComponents',NaN(4,1),'boundaryComponents',NaN(4,1), ...
    'witness',[],'degreeLossRadius',Inf,'degreeLossWitness',[], ...
    'message','','diagnostics',struct());
timer=tic;R.nominal=qstab.internal.isstableReport(parts(M.C0),'Variable','z','Grade',n,'Tolerance',o.Tolerance);
if ~strcmp(R.nominal.status,'stable')
    if strcmp(R.nominal.status,'undetermined'),R.status='undetermined';end
    R.message='The nominal polynomial must be numerically verified Schur at the declared forward grade.';return
end
if ~M.hasPerturbation
    R.status='unbounded';R.lower=Inf;R.upper=Inf;R.estimate=Inf;R.gap=0;
    R.message='The perturbation polynomial is identically zero.';
    R.diagnostics=struct('processedCells',0,'completeDomainCovered',true,'elapsedSeconds',toc(timer));return
end
upperBound=Inf;best=[];angle=NaN;
if any(M.S1(:,end)~=0)
    a=M.S0(:,end);b=M.S1(:,end);bs=max(abs(b));bb=b/bs;ib=[bb(1);-bb(2:4)]/(bb'*bb);
    if strcmp(side,'left'),ts=-mul(a,ib)/bs;else,ts=-mul(ib,a)/bs;end
    t=scalevec(ts,M.logScaleRatio);rd=norm(t);
    if isfinite(rd)&&rd>0
        [W,valid]=degreeWitness(M,t,side,o.Tolerance);
        if valid
            R.degreeLossRadius=rd;R.degreeLossWitness=W;
            upperBound=rd;best=W;
        end
    end
end
if n==0
    if isempty(best),R.status='undetermined';R.lower=0;R.upper=Inf;R.gap=Inf;return;end
    R=finish(R,upperBound,upperBound,best,angle,'converged','A nonzero constant loses regularity at the zero polynomial.');
    R.diagnostics=struct('processedCells',0,'completeDomainCovered',true,'elapsedSeconds',toc(timer));return
end
cells=[0,pi,0,0];processed=0;evaluations=0;
shortcut=struct('attempted',false,'accepted',false,'processedCells',0);
if o.MaxTime==0
    R=finish(R,0,upperBound,best,angle,'bounded','Zero time budget; the whole unit semicircle remains unresolved.');
    R.diagnostics=diagnostics(cells,processed,evaluations,shortcut,timer,o);return
end
seed=linspace(0,pi,max(5,2*n+1));
for j=1:numel(seed)
    if toc(timer)>=o.MaxTime,break;end
    E=qstab.internal.unitquotient(M,seed(j),side,o.Tolerance);evaluations=evaluations+1;
    [upperBound,best,angle]=update(E,upperBound,best,angle);
end
if isfinite(upperBound)
    target=max(0,upperBound-.8*(o.AbsTol+o.RelTol*upperBound));
    shortcut=qstab.internal.unitpositivity(M,target,side,min(1023,floor((o.MaxCells-1)/4)), ...
        o.MaxDepth,max(0,o.MaxTime-toc(timer)));
    processed=shortcut.processedCells;
    if shortcut.accepted
        lb=min(target,R.degreeLossRadius);
        R=finish(R,lb,upperBound,best,angle,'converged','A whole-unit-circle matrix Bernstein inequality supplies the numerical lower bound.');
        cells(4)=lb;R.diagnostics=diagnostics(cells,processed,evaluations,shortcut,timer,o);return
    end
end
[row,E]=makecell(M,0,pi,0,side,o.Tolerance);cells=row;processed=processed+1;evaluations=evaluations+1;
[upperBound,best,angle]=update(E,upperBound,best,angle);
state='bounded';reason='A cell, depth or time budget was reached.';
while true
    lowerBound=min(R.degreeLossRadius,min(cells(:,4)));
    if ~isfinite(lowerBound)&&~isfinite(upperBound)
        state='undetermined';reason='No representable finite witness or positive finite lower bound.';break
    end
    if lowerBound>upperBound
        if lowerBound-upperBound>1024*eps*max(upperBound,realmin)
            state='undetermined';reason='Numerical lower and upper bounds are inconsistent.';lowerBound=0;break
        end
        lowerBound=max(0,upperBound*(1-1024*eps));
    end
    if isfinite(upperBound)&&upperBound-lowerBound<=o.AbsTol+o.RelTol*upperBound
        state='converged';reason='Full-semicircle numerical bounds meet the requested tolerance.';break
    end
    if processed+2>o.MaxCells||toc(timer)>=o.MaxTime,break;end
    [~,idx]=min(cells(:,4));node=cells(idx,:);
    if node(3)>=o.MaxDepth,break;end
    % Frozen parent split; helpers are ordinary local functions, not closures.
    splitPoint=node(1)+(node(2)-node(1))/2;
    if splitPoint==node(1)||splitPoint==node(2)
        reason='Angular subdivision reached floating-point resolution.';break
    end
    [left,EL]=makecell(M,node(1),splitPoint,node(3)+1,side,o.Tolerance);
    [right,ER]=makecell(M,splitPoint,node(2),node(3)+1,side,o.Tolerance);
    cells(idx,:)=left;cells(end+1,:)=right; %#ok<AGROW>
    processed=processed+2;evaluations=evaluations+2;
    [upperBound,best,angle]=update(EL,upperBound,best,angle);
    [upperBound,best,angle]=update(ER,upperBound,best,angle);
end
R=finish(R,lowerBound,upperBound,best,angle,state,reason);
R.diagnostics=diagnostics(cells,processed,evaluations,shortcut,timer,o);
if ~R.diagnostics.partitionVerified
    R.status='undetermined';R.lower=0;R.gap=R.upper;
    R.message='Angular leaf coverage is inconsistent; no computed positive lower bound is returned.';
end
end
function [row,E]=makecell(M,l,r,depth,side,tol)
c=l+(r-l)/2;h=(r-l)/2;E=qstab.internal.unitquotient(M,c,side,tol);lb=0;
ks=0:M.grade;d0=sum(sqrt(sum(M.S0.^2,1)).*min(2,ks*h));d1=sum(sqrt(sum(M.S1.^2,1)).*min(2,ks*h));
m=E.scaledMinNominal;r0=E.scaledLower;
if E.valid&&m>d0&&r0>0
    denominator=m/r0+d1;
    if denominator==0,scaledBound=Inf;else,scaledBound=(m-d0)/denominator;end
    lb=scalevalue(max(0,scaledBound)*(1-512*eps*(M.grade+1)),M.logScaleRatio);
    if isnan(lb),lb=0;end
end
row=[l,r,depth,lb];
end
function [ub,best,angle]=update(E,ub,best,angle)
if E.valid&&isfinite(E.upper)&&~isempty(E.witness)&&E.upper<ub
    ub=E.upper;best=E.witness;angle=E.angle;
end
end
function [W,valid]=degreeWitness(M,t,side,tol)
Q=M.inputC0;D=zeros(size(Q));
for j=1:M.grade+1
    if strcmp(side,'left'),D(:,j)=mul(t,M.inputC1(:,j));else,D(:,j)=mul(M.inputC1(:,j),t);end
    Q(:,j)=Q(:,j)+D(:,j);
end
idx=M.grade+1;if strcmp(M.variable,'zinv'),idx=1;end
sc=norm(M.inputC0(:,idx))+norm(D(:,idx));res=norm(Q(:,idx))/max(sc,realmin);
valid=all(isfinite(Q(:)))&&all(isfinite(D(:)))&&res<=max(1024*eps,100*tol);
W=struct('kind','degreeLoss','valid',valid,'side',side,'variable',M.variable, ...
    'testVariable','z','grade',M.grade,'parameterComponents',t,'parameterNorm',norm(t), ...
    'coefficientComponents',Q,'perturbationComponents',D,'polynomial',parts(Q), ...
    'leadingIndexInInput',idx,'leadingResidual',res,'boundaryComponents',[], ...
    'forwardBoundaryComponents',[],'evaluationResidual',NaN,'validated',false);
end
function R=finish(R,lb,ub,W,angle,state,message)
R.status=state;R.lower=lb;R.upper=ub;R.estimate=ub;R.gap=ub-lb;R.witness=W;
R.frequency=angle;R.boundaryParameter=angle;
if isfinite(angle)
    if angle==0,R.boundaryArgument=1;elseif angle==pi,R.boundaryArgument=-1;else,R.boundaryArgument=exp(1i*angle);end
end
R.message=message;
if ~isempty(W)
    R.limitingEvent=W.kind;R.criticalParameterComponents=W.parameterComponents;
    if strcmp(W.kind,'finiteZero'),R.boundaryComponents=W.boundaryComponents;end
end
end
function D=diagnostics(cells,p,e,S,timer,o)
s=sortrows(cells,1);covered=~isempty(s)&&s(1,1)==0&&s(end,2)==pi&& ...
    all(s(:,2)>s(:,1))&&isequal(s(2:end,1),s(1:end-1,2));
D=struct('processedCells',p,'pointEvaluations',e,'leafCount',size(cells,1), ...
    'parameterDomain',[0,pi],'parameterMeaning','angle of forward z, radians', ...
    'cellColumns',{{'left','right','depth','lower'}},'cells',[], ...
    'partitionVerified',covered,'completeDomainCovered',covered, ...
    'bernsteinShortcut',S,'elapsedSeconds',toc(timer),'AbsTol',o.AbsTol,'RelTol',o.RelTol);
if o.StoreCells,D.cells=s;end
end
function tf=integer(x,l)
tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=l&&x==fix(x);
end
function y=scalevalue(x,l)
if x==0,y=0;elseif isinf(x),y=Inf;else,y=exp(log(x)+l);end
end
function y=scalevec(x,l)
y=zeros(size(x));for k=1:numel(x),y(k)=sign(x(k))*scalevalue(abs(x(k)),l);end
end
function z=mul(x,y)
z=[x(1)*y(1)-x(2:4)'*y(2:4);x(1)*y(2:4)+y(1)*x(2:4)+cross(x(2:4),y(2:4))];
end
function P=parts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
