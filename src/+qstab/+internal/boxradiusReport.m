function R = boxradiusReport(P,widths,varargin)
%BOXRADIUSREPORT Internal implementation returning the full report.
ip=inputParser;
addParameter(ip,'Variable','s',@localText);addParameter(ip,'Grade',[]);
addParameter(ip,'Method','zeroexclusion',@localText);
addParameter(ip,'AbsTol',1e-4,@localPositive);addParameter(ip,'RelTol',1e-3,@localPositive);
addParameter(ip,'Tolerance',1e-10,@(x)localPositive(x)&&x<1);
addParameter(ip,'InitialRadius',1,@localPositive);
addParameter(ip,'MaxIterations',32,@localInteger);
addParameter(ip,'MaxCells',8192,@localInteger);addParameter(ip,'MaxDepth',30,@localInteger);
addParameter(ip,'MaxTime',30,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&~isnan(x)&&x>=0);
addParameter(ip,'MaxRadius',1e6,@localPositive);
addParameter(ip,'BoundaryPoints',[]);
addParameter(ip,'MaxFacets',4096,@localInteger);
addParameter(ip,'StoreCertificates',false,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});o=ip.Results;
numericNames={'AbsTol','RelTol','Tolerance','InitialRadius','MaxIterations','MaxCells','MaxDepth','MaxTime','MaxRadius','MaxFacets'};
for jj=1:numel(numericNames)
    if isfield(o,numericNames{jj}),o.(numericNames{jj})=double(o.(numericNames{jj}));end
end
v=lower(char(o.Variable));method=lower(char(o.Method));
if ~any(strcmp(v,{'s','z','zinv'})),error('qstab:boxradius:Variable','Variable must be s, z or zinv.');end
if ~any(strcmp(method,{'zeroexclusion','lyapunov'})),error('qstab:boxradius:Method','Method must be zeroexclusion or lyapunov.');end
M=qstab.internal.boxmodel(P,widths,o.Grade);
if M.grade>64,error('qstab:boxradius:DegreeLimit','Current dense global box tests support grade at most 64.');end
if strcmp(v,'s'),seeds=zeros(4,1);else,seeds=[1 -1;zeros(3,2)];end
if ~isempty(o.BoundaryPoints)
    [a,b,c,d]=qstab.components(o.BoundaryPoints);q=[a(:)';b(:)';c(:)';d(:)'];
    if strcmp(v,'s'),on=all(q(1,:)==0);
    else,on=all(abs(sqrt(sum(q.^2,1))-1)<=128*eps);end
    if ~on,error('qstab:boxradius:Boundary','BoundaryPoints must lie on the imaginary axis/sphere or the unit sphere, as appropriate.');end
    seeds=[seeds q];
    o.BoundaryPoints=struct('a',a,'b',b,'c',c,'d',d);
end
R=struct('schema','qstab.boxradius.1','status','assumptionFailed', ...
    'method','nested-component-box-tests-with-zonotope-upper-witnesses', ...
    'variable',v,'grade',M.grade,'weights',M.weights, ...
    'lower',NaN,'upper',NaN,'estimate',NaN,'gap',NaN, ...
    'witness',[],'limitingEvent','none','degreeLossRadius',Inf, ...
    'degreeLossWitness',[],'nominal',[],'message','','options',o, ...
    'diagnostics',struct(),'validated',false, ...
    'assurance','floating-point whole-family decisions and primal/dual witnesses');
timer=tic;N=qstab.internal.isstableReport(localParts(M.C),'Variable',v,'Grade',M.grade,'Tolerance',o.Tolerance);R.nominal=N;
if ~strcmp(N.status,'stable')
    R.message='A strictly stable nominal polynomial at the declared grade is required.';return
end
if ~M.hasPerturbation
    R.status='unbounded';R.lower=Inf;R.upper=Inf;R.estimate=Inf;R.gap=0;
    R.message='All components are fixed; the nominal family is unchanged.';
    R.diagnostics=struct('iterations',0,'processedCells',0,'pointEvaluations',0,'elapsedSeconds',toc(timer));return
end
[degree,degreeW]=qstab.internal.boxdegree(M,v);R.degreeLossRadius=degree;R.degreeLossWitness=degreeW;
lo=0;hi=degree;best=degreeW;trials={};pointData={};processed=0;iterations=0;
lastStable=[];blocked=[];stop='iterationBudget';
if M.grade==0 && any(M.weights(:,1)==0 & M.C(:,1)~=0)
    R.status='unbounded';R.lower=Inf;R.upper=Inf;R.estimate=Inf;R.gap=0;
    R.message='A fixed nonzero constant component prevents the zero polynomial at every radius.';
    R.diagnostics=struct('iterations',0,'processedCells',0,'pointEvaluations',0,'elapsedSeconds',toc(timer));return
end
if M.grade==0&&isfinite(degree)
    R.lower=degree;R.upper=degree;R.estimate=degree;R.gap=0;R.status='converged';
    R.witness=degreeW;R.limitingEvent='degreeLoss';
    R.message='A nonzero constant has no finite zero; the first excluded member is the zero polynomial.';
    R.diagnostics=struct('iterations',0,'processedCells',0,'pointEvaluations',0,'elapsedSeconds',toc(timer));return
end
for j=1:size(seeds,2)
    if toc(timer)>=o.MaxTime,break;end
    E=qstab.internal.boxpoint(M,seeds(:,j),'right',o.Tolerance,o.MaxFacets,true);
    pointData{end+1}=E; %#ok<AGROW>
    if E.valid&&~isempty(E.witness)&&E.witness.valid&&E.upper<hi
        hi=E.upper;best=E.witness;
    end
end
searchRadius=min(o.InitialRadius,o.MaxRadius);
for iteration=1:o.MaxIterations
    if isfinite(hi)&&hi-lo<=o.AbsTol+o.RelTol*hi,stop='gapTolerance';break;end
    remaining=o.MaxCells-processed;
    if remaining<=0,stop='cellBudget';break;end
    if toc(timer)>=o.MaxTime,stop='timeBudget';break;end
    if isfinite(hi)
        cap=min(hi,o.MaxRadius);
        interior=sort(unique(blocked(blocked>lo & blocked<cap)));
        partition=[lo interior cap];[~,pick]=max(diff(partition));
        proposed=partition(pick)+(partition(pick+1)-partition(pick))/2;
        if cap<=lo,proposed=lo;end
    else
        proposed=max(searchRadius,lo*2);
        proposed=min(proposed,o.MaxRadius);
    end
    if ~(proposed>lo)&&isfinite(proposed),stop='explorationCap';break;end
    if ~isfinite(proposed),stop='unrepresentableScale';break;end
    T=qstab.internal.boxtrial(M,proposed,v,method,remaining,o.MaxDepth, ...
        max(0,o.MaxTime-toc(timer)),o.Tolerance,o.StoreCertificates);
    iterations=iterations+1;processed=processed+T.processedCells;
    record=struct('radius',proposed,'status',T.status,'processedCells',T.processedCells, ...
        'lowerBefore',lo,'upperBefore',hi,'witness',T.witness,'message',T.message, ...
        'result',T.result,'memberChecks',{T.memberChecks}, ...
        'lowerAfter',lo,'upperAfter',hi);
    if strcmp(T.status,'stable')
        lo=proposed;lastStable=T;
        if ~isfinite(hi),searchRadius=min(o.MaxRadius,2*proposed);end
    elseif ~isempty(T.witness)&&T.witness.valid
        candidate=T.witness.weightedNorm*(1+512*eps*(M.grade+1));
        if candidate<hi,hi=candidate;best=T.witness;end
        if candidate<=lo
            stop='inconsistentBounds';record.lowerAfter=lo;record.upperAfter=hi;
            trials{end+1}=record;break %#ok<AGROW>
        end
    else
        blocked(end+1)=proposed; %#ok<AGROW>
        if ~isfinite(hi)
            if proposed>=o.MaxRadius
                stop='unresolvedAtExplorationCap';record.lowerAfter=lo;record.upperAfter=hi;
                trials{end+1}=record;break %#ok<AGROW>
            end
            searchRadius=min(o.MaxRadius,2*proposed);
        end
    end
    record.lowerAfter=lo;record.upperAfter=hi;trials{end+1}=record; %#ok<AGROW>
end
R.lower=lo;R.upper=hi;R.estimate=hi;R.gap=hi-lo;R.witness=best;
if ~isempty(best),R.limitingEvent=best.kind;end
if strcmp(stop,'inconsistentBounds')
    R.status='undetermined';R.lower=0;R.upper=Inf;R.estimate=NaN;R.gap=Inf;
    R.message='Contradictory numerical evidence: no global bracket is asserted.';
elseif isfinite(hi)&&hi-lo<=o.AbsTol+o.RelTol*hi
    R.status='converged';stop='gapTolerance';R.message='Numerical global bracket meets the requested gap tolerance.';
elseif lo>0||isfinite(hi)
    R.status='bounded';R.message=['Conservative global bracket retained; stop: ',stop,'.'];
else
    R.status='undetermined';R.message=['No positive global lower or finite upper found; stop: ',stop,'.'];
end
R.diagnostics=struct('iterations',iterations,'processedCells',processed, ...
    'pointEvaluations',numel(pointData),'pointData',{pointData},'trials',{trials}, ...
    'lastStableCertificate',lastStable,'stopReason',stop,'elapsedSeconds',toc(timer), ...
    'globalLowerFromWholeFamilies',true,'directDirectionMinimumClaim',false, ...
    'explorationCap',o.MaxRadius,'unresolvedTrialRadii',blocked, ...
    'unresolvedDoesNotUpdateBounds',true);
end
function P=localParts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
function y=localText(x),y=ischar(x)||(isstring(x)&&isscalar(x));end
function y=localPositive(x),y=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0;end
function y=localInteger(x),y=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x==fix(x);end
