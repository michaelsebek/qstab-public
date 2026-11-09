function D = certificateDataReport(R)
%CERTIFICATEDATAREPORT Internal implementation returning the full report.
if ~isstruct(R)||~isscalar(R)||~all(isfield(R,{'status','method','variable','certificates','unresolvedRegions'}))
    error('qstab:certificateData:Report','Supply a scalar qstab robust-test report with retained-cell fields.');
end
if ~iscell(R.certificates)||~iscell(R.unresolvedRegions)|| ...
        ~ismember(R.status,{'stable','notStable','undetermined','assumptionFailed'})
    error('qstab:certificateData:Report','Invalid report status or cell containers.');
end
tv=R.variable;if strcmp(tv,'zinv'),tv='z';end
if isfield(R,'testVariable'),tv=R.testVariable;end
D=struct('schema','qstab.certificate-data.1','status',R.status,'method',R.method, ...
    'variable',R.variable,'testVariable',tv,'reportedComplete',false, ...
    'storageRequested',false,'storageComplete',false,'storedAcceptedCount',numel(R.certificates), ...
    'unresolvedDescriptorCount',numel(R.unresolvedRegions),'reportedAcceptedCount',NaN, ...
    'reportedUnresolvedCount',NaN,'hasWitness',false,'witness',[], ...
    'rows',repmat(localEmpty(),0,1),'reduction',[], ...
    'isCertificate',false,'validated',false,'description','');
if isfield(R,'certificateComplete'),D.reportedComplete=isequal(R.certificateComplete,true);end
if isfield(R,'options')&&isfield(R.options,'StoreCertificates'),D.storageRequested=isequal(R.options.StoreCertificates,true);end
if isfield(R,'certifiedCells'),D.reportedAcceptedCount=R.certifiedCells;end
if isfield(R,'unresolvedCells'),D.reportedUnresolvedCount=R.unresolvedCells;end
D.storageComplete=D.storageRequested&&isequal(D.reportedAcceptedCount,D.storedAcceptedCount);
if isfield(R,'witness'),D.witness=R.witness;D.hasWitness=~isempty(R.witness);end
if isfield(R,'reduction'),D.reduction=R.reduction;end
for j=1:numel(R.certificates)
    c=R.certificates{j};
    if ~isstruct(c)||~isscalar(c)||~isfield(c,'cell')||~isstruct(c.cell)||~isscalar(c.cell)
        error('qstab:certificateData:Certificate','Stored certificate %d is malformed.',j);
    end
    row=localRead(c.cell);row.state='accepted';row.sourceIndex=j;row.record=c;
    if isfield(c,'generatorIndices'),row.generatorIndices=c.generatorIndices;end
    if isfield(c,'faceIndex'),row.faceIndex=c.faceIndex;end
    if isfield(c,'slack')
        row.metric=c.slack;row.metricName='Lyapunov inequality slack';
        row.metricUnits='normalizedLyapunovInequality';row.reason=localField(c,'reason','');
    elseif isfield(c,'diagnostic')&&isstruct(c.diagnostic)
        d=c.diagnostic;row.reason=localField(d,'reason',row.reason);
        if isfield(d,'gap')
            row.metric=d.gap;row.metricName='Stored norm-value separation gap';
        elseif all(isfield(d,{'centerDistance','freqVariation'}))
            row.metric=d.centerDistance-d.freqVariation;
            row.metricName='Distance minus frequency variation (before tolerance allowance)';
        elseif all(isfield(d,{'realInterval','imagInterval'}))
            x=d.realInterval;y=d.imagInterval;
            row.metric=max([x(1),-x(2),y(1),-y(2)]);
            row.metricName='Axis separation of interval enclosure (before tolerance allowance)';
        end
        row.metricUnits='normPolynomialValue';
        row.enclosureAvailable=all(isfield(d,{'realRange','imagRange'}))|| ...
            all(isfield(d,{'realInterval','imagInterval'}))|| ...
            (~isempty(row.vertexWeights)&&all(isfield(d,{'centerDistance','freqVariation'})));
    end
    D.rows(end+1,1)=row; %#ok<AGROW>
end
for j=1:numel(R.unresolvedRegions)
    r=R.unresolvedRegions{j};
    if ~isstruct(r)||~isscalar(r)
        error('qstab:certificateData:Region','Unresolved descriptor %d is malformed.',j);
    end
    row=localRead(r);row.state='unresolved';row.sourceIndex=j;row.record=r;
    D.rows(end+1,1)=row; %#ok<AGROW>
end
D.description=['Read-only diagnostic of stored records. Frequency bars are projections, ', ...
    'not an independent proof of coverage. Missing stored leaves, unvisited ', ...
    'subfamilies and a reported solver verdict remain distinct.'];
end
function r=localEmpty()
r=struct('state','','sourceIndex',0,'kind','','depth',NaN,'reason','', ...
    'frequencyBounds',[],'parameterBounds',[],'vertexWeights',[], ...
    'generatorIndices',[],'faceIndex',[],'evenBox',[],'oddBox',[], ...
    'metric',NaN,'metricName','','metricUnits','none', ...
    'enclosureAvailable',false,'record',[]);
end
function r=localRead(x)
r=localEmpty();r.kind=localField(x,'kind',localField(x,'uncertainty','parameterCell'));
r.depth=localField(x,'depth',NaN);r.reason=localField(x,'reason','');
fields={'parameterBounds','vertexWeights','generatorIndices','faceIndex','evenBox','oddBox'};
for j=1:numel(fields),r.(fields{j})=localField(x,fields{j},[]);end
if isempty(r.parameterBounds)&&all(isfield(x,{'lower','upper'})),r.parameterBounds=[x.lower,x.upper];end
if isempty(r.vertexWeights)&&isfield(x,'weights'),r.vertexWeights=x.weights;end
f=localField(x,'frequencyBounds',[]);
if isnumeric(f)&&isreal(f)&&numel(f)==2&&all(isfinite(f))&&f(1)<=f(2),r.frequencyBounds=reshape(f,1,2);end
end
function x=localField(s,f,d)
if isfield(s,f),x=s.(f);else,x=d;end
end
