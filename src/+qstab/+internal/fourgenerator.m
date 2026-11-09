function R = fourgenerator(F,o,variable)
%FOURGENERATOR Internal qB03 reduction; called only by qstab.zeroexclude.
%   Enumerates C(m,4) generator subfamilies lazily, or the complete family
%   if m<=4. Each subfamily is tested by the unchanged Bernstein zero-
%   exclusion engine. It is NOT a vertex test and does not construct the
%   exposed faces of a coefficient polytope.
%
%   qB03 Corollary "Four-generator-subfamily test" and Theorem "Schur
%   three-face counterpart" require invariant degree over the FULL family.
%   The nonzero distance of the FULL leading hull is certified first.
%   Testing four-generator leading hulls alone would not rule out degree
%   loss inside a five-generator hull. No coefficient normalization is used.
%
%   No normalization, independent coefficient enclosure, new class, or
%   external optimizer is used. All output weights refer to original input
%   generator indices, including duplicates and zero high-degree padding.

started=tic;
testVariable=variable;
if strcmp(variable,'zinv'), testVariable='z'; end
R=struct('status','undetermined','isStable',[], ...
    'method','fourGeneratorZeroExclusion','uncertainty',F.kind, ...
    'variable',variable,'testVariable',testVariable,'degree',F.degree, ...
    'message','','assurance','floatingPointUndetermined', ...
    'validated',false,'certificateComplete',false,'witness',[], ...
    'nominal',[],'frequencyLimit',[],'processedCells',0,'certifiedCells',0, ...
    'unresolvedCells',0,'maximumDepth',0,'elapsedTime',0, ...
    'unresolvedRegions',{{}},'certificates',{{}},'options',o);
R.reduction=struct('name','four-generator','fullVertexCount',0, ...
    'subsetSize',0,'parameterDimensionBound',0,'totalSubfamilies',0, ...
    'totalCountExact',true,'processedSubfamilies',0,'stableSubfamilies',0, ...
    'unresolvedSubfamilies',0,'unvisitedSubfamilies',0, ...
    'nextGeneratorIndices',[],'allVisited',false,'complete',false, ...
    'enumeration','lexicographic, generated lazily','subfamilies',{{}}, ...
    'referenceGeneratorIndices',[], ...
    'degreeAssumption','full leading-coefficient hull excludes zero in tested variable', ...
    'theorem','qB03 four-generator subfamilies; Schur counterpart after delay reversal');
if ~strcmp(F.kind,'polytope')
    R.message='Four-generator reduction supports scalar polynomial polytopes only.';
    R.unresolvedRegions={struct('kind','unsupportedReduction', ...
        'originalUncertainty',F.kind,'reason',R.message)};
    R=localFinish(R,started,[]);
    return
end

V=F.vertices; m=size(V,3); k=min(m,4);
[total,countExact]=localCount(m,k);
R.reduction.fullVertexCount=m;
R.reduction.subsetSize=k;
R.reduction.parameterDimensionBound=k-1;
R.reduction.totalSubfamilies=total;
R.reduction.totalCountExact=countExact;
next=1:k;
leadIndex=F.degree+1;
if strcmp(variable,'zinv'), leadIndex=1; end
for j=1:m
    if all(V(:,leadIndex,j)==0)
        weights=zeros(m,1); weights(j)=1;
        C=V(:,:,j); Ct=C;
        if strcmp(variable,'zinv'), Ct=C(:,end:-1:1); end
        R.status='assumptionFailed';
        R.assurance='floatingPointAssumptionWitness';
        R.witness=struct('components',Ct,'originalComponents',C, ...
            'originalWeights',weights,'kind','zeroLeadingCoefficient', ...
            'coefficientVariable',testVariable, ...
            'originalCoefficientVariable',variable,'generatorIndices',j);
        R.message='An admissible generator has zero leading coefficient in the tested variable.';
        R=localPending(R,next,R.message);
        R=localFinish(R,started,next);
        return
    end
end
R.leading=qstab.internal.leadbound(F,variable,o.Tolerance,256, ...
    max(0,o.MaxTime-toc(started)));
if strcmp(R.leading.status,'zeroWitness')
    R.status='assumptionFailed';R.assurance='floatingPointAssumptionWitness';
    R.witness=R.leading.witness;R.message=R.leading.message;
    R=localPending(R,next,R.message);R=localFinish(R,started,next);return
elseif ~R.leading.isSeparated
    R.message=['Full leading hull is unresolved before reduction: ',R.leading.message];
    R=localPending(R,next,R.message);R=localFinish(R,started,next);return
end

reason='';
while ~isempty(next)
    if R.reduction.processedSubfamilies>=o.MaxSubfamilies
        reason='MaxSubfamilies reached'; break
    end
    remainingCells=o.MaxCells-R.processedCells;
    if remainingCells<1
        reason='Global MaxCells reached'; break
    end
    remainingTime=o.MaxTime-toc(started);
    if remainingTime<=0
        reason='Global MaxTime reached'; break
    end
    indices=next;
    % Rebuild from PADDED input components so all subfamilies retain the
    % original grade, notably when some delay generators end with zeros.
    data=cell(1,k);
    for j=1:k, data{j}=localPolynomial(V(:,:,indices(j))); end
    remainingTime=max(0,o.MaxTime-toc(started));
    if remainingTime==0
        reason='Global MaxTime reached'; break
    end
    child=qstab.internal.zeroexcludeReport(data,'Uncertainty','polytope', ...
        'Variable',variable,'Reduction','none', ...
        'MaxCells',remainingCells,'MaxDepth',o.MaxDepth, ...
        'MaxTime',remainingTime,'Tolerance',o.Tolerance, ...
        'StoreCertificates',o.StoreCertificates);

    R.reduction.processedSubfamilies=R.reduction.processedSubfamilies+1;
    R.processedCells=R.processedCells+child.processedCells;
    R.certifiedCells=R.certifiedCells+child.certifiedCells;
    R.maximumDepth=max(R.maximumDepth,child.maximumDepth);
    if ~isempty(child.frequencyLimit)
        if isempty(R.frequencyLimit), R.frequencyLimit=child.frequencyLimit;
        else, R.frequencyLimit=max(R.frequencyLimit,child.frequencyLimit); end
    end
    if isempty(R.nominal) && ~isempty(child.nominal)
        R.nominal=child.nominal;
        R.reduction.referenceGeneratorIndices=indices;
    end
    entry=struct('generatorIndices',indices,'status',child.status, ...
        'certificateComplete',child.certificateComplete, ...
        'processedCells',child.processedCells,'certifiedCells',child.certifiedCells, ...
        'unresolvedCells',child.unresolvedCells,'maximumDepth',child.maximumDepth, ...
        'frequencyLimit',child.frequencyLimit,'message',child.message);
    R.reduction.subfamilies{end+1,1}=entry;
    for j=1:numel(child.certificates)
        cert=child.certificates{j};
        cert.cell=localMapRegion(cert.cell,indices,m);
        cert.generatorIndices=indices;
        R.certificates{end+1,1}=cert; %#ok<AGROW>
    end
    for j=1:numel(child.unresolvedRegions)
        R.unresolvedRegions{end+1,1}=localMapRegion( ...
            child.unresolvedRegions{j},indices,m); %#ok<AGROW>
    end
    next=localNext(indices,m);

    if strcmp(child.status,'stable') && child.certificateComplete
        R.reduction.stableSubfamilies=R.reduction.stableSubfamilies+1;
    elseif strcmp(child.status,'notStable') || strcmp(child.status,'assumptionFailed')
        R.status=child.status; R.isStable=child.isStable;
        R.assurance=child.assurance;
        R.witness=localMapWitness(child.witness,indices,V,variable,testVariable);
        R.message=['An admissible original-family witness was found in a ', ...
            'four-generator subfamily: ',child.message];
        R=localPending(R,next,'Not visited after an admissible witness');
        R=localFinish(R,started,next);
        return
    else
        R.reduction.unresolvedSubfamilies=R.reduction.unresolvedSubfamilies+1;
        % Retain ambiguity even if another subfamily subsequently succeeds.
        if isempty(child.unresolvedRegions)
            region=struct('kind','unresolvedGeneratorSubfamily', ...
                'generatorIndices',indices,'fullVertexCount',m, ...
                'vertexWeights',localMapWeights(eye(k),indices,m), ...
                'reason',child.message);
            R.unresolvedRegions{end+1,1}=region; %#ok<AGROW>
        end
    end
end

if isempty(next) && R.reduction.unresolvedSubfamilies==0 && ...
        R.reduction.stableSubfamilies==R.reduction.processedSubfamilies
    R.status='stable'; R.isStable=true; R.certificateComplete=true;
    R.assurance='floatingPointSufficientCertificate';
    R.message=['Every required complete generator subfamily is stable by ', ...
        'norm-value zero exclusion; the qB03 four-generator theorem applies.'];
    R.reduction.complete=true;
else
    R=localPending(R,next,reason);
    if isempty(reason), reason='One or more generator subfamilies remain unresolved'; end
    R.message=['No complete four-generator certificate: ',reason];
end
R=localFinish(R,started,next);
end

function S=localPolynomial(C)
S=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end

function W=localMapWeights(w,indices,m)
W=zeros(m,size(w,2));
W(indices,:)=w;
end

function region=localMapRegion(region,indices,m)
if isfield(region,'vertexWeights') && ~isempty(region.vertexWeights)
    region.vertexWeights=localMapWeights(region.vertexWeights,indices,m);
end
region.generatorIndices=indices;
region.fullVertexCount=m;
end

function W=localMapWitness(W,indices,V,variable,testVariable)
if isempty(W), return; end
W.generatorIndices=indices;
if isfield(W,'originalWeights') && ~isempty(W.originalWeights)
    W.subfamilyWeights=W.originalWeights;
    W.originalWeights=localMapWeights(W.originalWeights,indices,size(V,3));
    original=zeros(size(V,1),size(V,2));
    for j=1:numel(indices)
        original=original+W.subfamilyWeights(j)*V(:,:,indices(j));
    end
    W.originalComponents=original;
end
% Existing child.components remain in the variable actually tested. The
% separate originalComponents remove ambiguity for delay witnesses.
W.coefficientVariable=testVariable;
W.originalCoefficientVariable=variable;
end

function R=localPending(R,next,reason)
if isempty(next), return; end
if R.reduction.totalCountExact
    count=R.reduction.totalSubfamilies-R.reduction.processedSubfamilies;
else
    count=Inf;
end
% A compact exact descriptor of an entire lexicographic suffix. Avoid
% allocating nchoosek(1:m,4), which can be enormous even for a small budget.
region=struct('kind','unvisitedGeneratorSubfamilies', ...
    'nextGeneratorIndices',next,'subsetSize',R.reduction.subsetSize, ...
    'fullVertexCount',R.reduction.fullVertexCount, ...
    'count',count,'countExact',R.reduction.totalCountExact, ...
    'enumeration','lexicographic suffix including nextGeneratorIndices', ...
    'reason',reason);
R.unresolvedRegions{end+1,1}=region;
end

function R=localFinish(R,started,next)
R.reduction.nextGeneratorIndices=next;
R.reduction.allVisited=isempty(next) && R.reduction.processedSubfamilies>0;
if isempty(next)
    R.reduction.unvisitedSubfamilies=0;
elseif R.reduction.totalCountExact
    R.reduction.unvisitedSubfamilies= ...
        R.reduction.totalSubfamilies-R.reduction.processedSubfamilies;
else
    R.reduction.unvisitedSubfamilies=Inf;
end
% One compact suffix descriptor may stand for many unvisited subfamilies.
% unresolvedCells counts stored region descriptors, not unallocated cells.
R.unresolvedCells=numel(R.unresolvedRegions);
R.elapsedTime=toc(started);
end

function next=localNext(indices,m)
next=indices; k=numel(indices); j=k;
while j>=1 && next(j)==m-k+j, j=j-1; end
if j==0, next=[]; return; end
next(j)=next(j)+1;
for ell=j+1:k, next(ell)=next(ell-1)+1; end
end

function [count,isExact]=localCount(m,k)
% Cancel the small factorial before multiplication. Products remain exact
% integer doubles while the final count is <= flintmax. The count is only
% diagnostic: loop termination is always determined by localNext.
factors=m-k+(1:k);
for denom=2:k
    left=denom;
    for j=1:k
        common=gcd(factors(j),left);
        factors(j)=factors(j)/common; left=left/common;
        if left==1, break; end
    end
end
count=prod(factors);
isExact=isfinite(count) && count<=flintmax;
if ~isExact, count=Inf; end
end
