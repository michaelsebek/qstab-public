function R = intervalreduce(F,o,variable,mode)
%INTERVALREDUCE qB02 cube, pair-free faces and monic-quadratic two-faces.
% MaxCells/MaxTime/MaxSubfamilies are GLOBAL budgets. Every face is a whole
% continuous parameter family, never just a collection of endpoints.
started=tic;
R=struct('status','undetermined','isStable',[], ...
    'method','intervalFaceZeroExclusion','uncertainty',F.kind, ...
    'variable',variable,'testVariable',variable,'degree',F.degree, ...
    'message','','assurance','floatingPointUndetermined','validated',false, ...
    'certificateComplete',false,'witness',[],'nominal',[], ...
    'frequencyLimit',[],'processedCells',0,'certifiedCells',0, ...
    'unresolvedCells',0,'maximumDepth',0,'elapsedTime',0, ...
    'unresolvedRegions',{{}},'certificates',{{}},'options',o);
R.reduction=struct('name',mode,'parameterDimensionBound',0, ...
    'totalSubfamilies',0,'processedSubfamilies',0,'stableSubfamilies',0, ...
    'unresolvedSubfamilies',0,'unvisitedSubfamilies',0,'nextFaceIndex',[], ...
    'complete',false,'allVisited',false,'subfamilies',{{}}, ...
    'theorem','','vertexOnlyTest',false);
if ~strcmp(F.kind,'interval') || ~strcmp(variable,'s')
    R.message='The selected qB02 reduction requires a scalar Cartesian interval in s.';
    R.unresolvedRegions={struct('kind','unsupportedReduction','reason',R.message)};
    R=localFinish(R,started,1);return
end
if strcmp(mode,'quadratic-faces')
    unit=[1;0;0;0];
    if F.degree~=2 || ~isequal(F.intervalLower(:,end),unit) || ...
            ~isequal(F.intervalUpper(:,end),unit)
        R.message='quadratic-faces requires exactly monic s^2+a*s+b; no family normalization is performed.';
        R.unresolvedRegions={struct('kind','unsupportedReduction','reason',R.message)};
        R=localFinish(R,started,1);return
    end
end
R.leading=qstab.internal.leadbound(F,variable,o.Tolerance,256, ...
    max(0,o.MaxTime-toc(started)));
if strcmp(R.leading.status,'zeroWitness')
    R.status='assumptionFailed';R.assurance='floatingPointAssumptionWitness';
    R.witness=R.leading.witness;R.message=R.leading.message;
    R=localFinish(R,started,1);return
elseif ~R.leading.isSeparated
    R.message=['Full interval degree condition is unresolved: ',R.leading.message];
    R.unresolvedRegions={struct('kind','leadingCoefficientUnresolved','reason',R.message)};
    R=localFinish(R,started,1);return
end
if strcmp(mode,'quadratic-faces')
    A=F;A.kind='affine';tasks=localQuadraticTasks(F);
    R.reduction.parameterDimensionBound=min(2,sum(F.parameterIndex(:,2)==1));
    R.reduction.theorem='qB02 monic quadratic interval theorem: two-faces in a and minimal-real vertices in b';
else
    [A,cubeInfo]=qstab.internal.cubefamily(F);R.reduction.cube=cubeInfo;
    if ~cubeInfo.representable
        R.message='Kharitonov endpoint differences exceed the representable range.';
        R.unresolvedRegions={struct('kind','unrepresentableCube','reason',R.message)};
        R=localFinish(R,started,1);return
    end
    if strcmp(mode,'kharitonov-cube')
        tasks={struct('bounds',A.bounds,'freeParameters',cubeInfo.activeParameters, ...
            'fixedParameters',setdiff(1:8,cubeInfo.activeParameters))};
        R.reduction.parameterDimensionBound=cubeInfo.activeParameterCount;
        R.reduction.theorem='qB02 exact eight-parameter Kharitonov cube';
    else
        tasks=localPairFreeTasks(A,cubeInfo);
        R.reduction.parameterDimensionBound=min(3,numel(unique( ...
            cubeInfo.parameterComponents(cubeInfo.activeParameters))));
        R.reduction.theorem='qB02 pair-free faces of the Kharitonov cube (at most 1024)';
    end
end
R.reduction.totalSubfamilies=numel(tasks);next=1;reason='';
while next<=numel(tasks)
    if R.reduction.processedSubfamilies>=o.MaxSubfamilies
        reason='Global MaxSubfamilies reached';break
    end
    cells=o.MaxCells-R.processedCells;remaining=o.MaxTime-toc(started);
    if cells<1,reason='Global MaxCells reached';break;end
    if remaining<=0,reason='Global MaxTime reached';break;end
    T=tasks{next};childF=A;childF.bounds=T.bounds;
    childOpts=o;childOpts.Reduction='none';childOpts.MaxCells=cells;childOpts.MaxTime=remaining;
    child=qstab.internal.affinezeroexclude(childF,childOpts,variable);
    R.reduction.processedSubfamilies=R.reduction.processedSubfamilies+1;
    R.processedCells=R.processedCells+child.processedCells;
    R.certifiedCells=R.certifiedCells+child.certifiedCells;
    R.maximumDepth=max(R.maximumDepth,child.maximumDepth);
    if ~isempty(child.frequencyLimit)
        if isempty(R.frequencyLimit),R.frequencyLimit=child.frequencyLimit;
        else,R.frequencyLimit=max(R.frequencyLimit,child.frequencyLimit);end
    end
    if isempty(R.nominal),R.nominal=child.nominal;end
    entry=struct('faceIndex',next,'freeParameters',T.freeParameters, ...
        'fixedParameters',T.fixedParameters,'parameterBounds',T.bounds, ...
        'status',child.status,'certificateComplete',child.certificateComplete, ...
        'processedCells',child.processedCells,'certifiedCells',child.certifiedCells, ...
        'unresolvedCells',child.unresolvedCells,'message',child.message);
    R.reduction.subfamilies{end+1,1}=entry;
    for j=1:numel(child.certificates)
        c=child.certificates{j};c.faceIndex=next;c.reduction=mode;
        R.certificates{end+1,1}=c; %#ok<AGROW>
    end
    for j=1:numel(child.unresolvedRegions)
        c=child.unresolvedRegions{j};c.faceIndex=next;c.reduction=mode;
        R.unresolvedRegions{end+1,1}=c; %#ok<AGROW>
    end
    current=next;next=next+1;
    if strcmp(child.status,'stable') && child.certificateComplete
        R.reduction.stableSubfamilies=R.reduction.stableSubfamilies+1;
    elseif strcmp(child.status,'notStable')
        W=child.witness;W.faceIndex=current;W.reduction=mode;
        if ~strcmp(mode,'quadratic-faces'),W.cubeParameter=W.parameter;end
        % Coefficients are an actual member, not a Bernstein control point.
        C=W.originalComponents;L=F.intervalLower;U=F.intervalUpper;
        slack=o.Tolerance*max(1,max(abs([L(:);U(:)])));
        if all(C(:)>=L(:)-slack) && all(C(:)<=U(:)+slack)
            R.status='notStable';R.isStable=false;R.assurance='floatingPointCounterexample';
            R.witness=W;R.message='An actual original-interval member in a required face is not stable.';
            R=localPending(R,tasks,next,'Unvisited after an actual interval witness');
            R=localFinish(R,started,next);return
        end
        R.reduction.unresolvedSubfamilies=R.reduction.unresolvedSubfamilies+1;
        R.unresolvedRegions{end+1,1}=struct('kind','witnessReconstructionUnresolved', ...
            'faceIndex',current,'parameterBounds',T.bounds, ...
            'reason','Candidate could not be verified in the original coefficient interval');
    else
        R.reduction.unresolvedSubfamilies=R.reduction.unresolvedSubfamilies+1;
        if isempty(child.unresolvedRegions)
            R.unresolvedRegions{end+1,1}=struct('kind','unresolvedIntervalFace', ...
                'faceIndex',current,'parameterBounds',T.bounds,'reason',child.message); %#ok<AGROW>
        end
    end
end
if next>numel(tasks) && R.reduction.unresolvedSubfamilies==0
    R.status='stable';R.isStable=true;R.certificateComplete=true;
    R.assurance='floatingPointSufficientCertificate';R.reduction.complete=true;
    R.message='Every required continuous qB02 subfamily has a complete norm-value certificate.';
else
    if isempty(reason),reason='One or more required interval subfamilies remain unresolved';end
    R=localPending(R,tasks,next,reason);R.message=['Incomplete qB02 reduction: ',reason];
end
R=localFinish(R,started,next);
end

function tasks=localPairFreeTasks(A,info)
active=info.activeParameters;labels=info.parameterComponents;
components=unique(labels(active));k=min(3,numel(components));
sets=localCombinations(components,k);tasks={};
for row=1:size(sets,1)
    chosen=sets(row,:);choices={zeros(1,0)};
    for nu=chosen
        cand=active(labels(active)==nu);out={};
        for j=1:numel(choices)
            for p=cand,out{end+1}=[choices{j},p];end %#ok<AGROW>
        end
        choices=out;
    end
    for j=1:numel(choices)
        free=sort(choices{j});fixed=setdiff(active,free);
        for mask=0:2^numel(fixed)-1
            B=zeros(8,2);B(free,2)=1;
            for ell=1:numel(fixed)
                bit=bitget(mask,ell);B(fixed(ell),:)=[bit,bit];
            end
            tasks{end+1}=struct('bounds',B,'freeParameters',free, ...
                'fixedParameters',setdiff(1:8,free)); %#ok<AGROW>
        end
    end
end
if isempty(active)
    tasks={struct('bounds',zeros(size(A.bounds)),'freeParameters',[], ...
        'fixedParameters',1:8)};
end
end

function tasks=localQuadraticTasks(F)
a=find(F.parameterIndex(:,2)==1).';
bIm=find(F.parameterIndex(:,2)==0 & F.parameterIndex(:,1)>1).';
bReal=find(F.parameterIndex(:,2)==0 & F.parameterIndex(:,1)==1).';
k=min(2,numel(a));sets=localCombinations(a,k);tasks={};
for row=1:size(sets,1)
    free=sets(row,:);aFixed=setdiff(a,free);fixed=[aFixed,bIm];
    for mask=0:2^numel(fixed)-1
        B=F.bounds;
        for idx=bReal,B(idx,2)=B(idx,1);end
        for ell=1:numel(fixed)
            idx=fixed(ell);v=F.bounds(idx,1+bitget(mask,ell));B(idx,:)=[v,v];
        end
        tasks{end+1}=struct('bounds',B,'freeParameters',free, ...
            'fixedParameters',setdiff(1:size(B,1),free)); %#ok<AGROW>
    end
end
end
function C=localCombinations(v,k)
if k==0,C=zeros(1,0);else,C=nchoosek(v,k);end
end
function R=localPending(R,tasks,next,reason)
if next>numel(tasks),return;end
% At most 1024 descriptors; retain exact remaining bounds, not sampled faces.
R.unresolvedRegions{end+1,1}=struct('kind','unvisitedIntervalFaces', ...
    'nextFaceIndex',next,'count',numel(tasks)-next+1, ...
    'faces',{tasks(next:end)},'reason',reason);
end
function R=localFinish(R,started,next)
R.reduction.unvisitedSubfamilies=max(0,R.reduction.totalSubfamilies-R.reduction.processedSubfamilies);
if next<=R.reduction.totalSubfamilies,R.reduction.nextFaceIndex=next;
else,R.reduction.nextFaceIndex=[];end
R.reduction.allVisited=R.reduction.totalSubfamilies>0 && R.reduction.unvisitedSubfamilies==0;
R.unresolvedCells=numel(R.unresolvedRegions);R.elapsedTime=toc(started);
end
