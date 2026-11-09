function E = certificateEnclosureReport(data,R,varargin)
%CERTIFICATEENCLOSUREREPORT Internal implementation returning the full report.
ip=inputParser;ip.FunctionName='qstab.certificateEnclosure';
addParameter(ip,'CellIndex',1,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==fix(x));
addParameter(ip,'DiskSides',64,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=8&&x<=2048&&x==fix(x));
parse(ip,varargin{:});o=ip.Results;D=qstab.internal.certificateDataReport(R);
if o.CellIndex>D.storedAcceptedCount,error('qstab:certificateEnclosure:Index','The requested accepted leaf is not stored.');end
row=D.rows(o.CellIndex);c=row.record;
if ~row.enclosureAvailable||~isfield(c,'diagnostic')
    error('qstab:certificateEnclosure:NoGeometry','This leaf has no stored norm-value enclosure; use certificateData/plotCertificate.');
end
d=c.diagnostic;f=row.frequencyBounds;
if isempty(f),error('qstab:certificateEnclosure:Frequency','The leaf has no finite frequency interval.');end
mid=.5*f(1)+.5*f(2);if strcmp(D.testVariable,'s'),z=1i*mid;else,z=exp(1i*mid);end
radius=0;controls=[];binding='geometryReadFromReport';residual=NaN;
if ~isempty(row.vertexWeights)
    if isempty(data),error('qstab:certificateEnclosure:Family','The original polynomial polytope is required.');end
    F=qstab.familycomponents(data,'Uncertainty','polytope');V=F.vertices;
    if ~isequal(F.degree,R.degree),error('qstab:certificateEnclosure:FamilyMismatch','Original grade differs from the report.');end
    if strcmp(D.variable,'zinv'),V=V(:,end:-1:1,:);end
    polynomials=cell(1,size(V,3));
    for j=1:numel(polynomials),polynomials{j}=localPolynomial(V(:,:,j));end
    base=qstab.internal.valueenclosureReport(polynomials,z,'VertexCells',{row.vertexWeights});
    p=base.regions{1}.polygon;controls=base.regions{1}.controls;
    % Check reconstructed distance and derivative variation against storage.
    m=size(V,3);W=base.regions{1}.vertexWeights;
    C=reshape(reshape(V,[],m)*W,4,F.degree+1,size(W,2));n=F.degree;L=0;
    if strcmp(D.testVariable,'s'),rho=max(1,f(2));else,rho=1;end
    for a=1:size(C,3)
        for b=a:size(C,3)
            beta=zeros(1,2*n+1);
            for nu=1:4,beta=beta+conv(C(nu,:,a),C(nu,:,b));end
            L=max(L,sum((1:2*n).*abs(beta(2:end)).*rho.^(0:2*n-1)));
        end
    end
    computedRadius=L*(f(2)-f(1))/2;
    computedDistance=qstab.internal.diagnosticGeometry('distance',p);
    residual=max(abs([computedDistance-d.centerDistance,computedRadius-d.freqVariation]));
    scale=max([1;abs(controls(:));abs(computedRadius);abs(d.freqVariation)]);
    tol=1e-10;if isfield(R,'options')&&isfield(R.options,'Tolerance'),tol=R.options.Tolerance;end
    if ~isfinite(residual)||residual>max(1024*eps,10*tol)*scale
        error('qstab:certificateEnclosure:FamilyMismatch','Reconstructed geometry disagrees with the saved leaf.');
    end
    radius=d.freqVariation;binding='originalWeightsAndSavedGeometryChecked';kind='BernsteinPolygonPlusFrequencyDisk';
elseif all(isfield(d,{'realRange','imagRange'}))
    p=qstab.internal.diagnosticGeometry('rectangle',[d.realRange;d.imagRange]);
    radius=d.freqVariation;kind='JointRectanglePlusFrequencyDisk';
elseif all(isfield(d,{'realInterval','imagInterval'}))
    p=qstab.internal.diagnosticGeometry('rectangle',[d.realInterval;d.imagInterval]);
    kind='WholeFrequencyEvenOddRectangle';
else
    error('qstab:certificateEnclosure:NoGeometry','Unsupported stored enclosure fields.');
end
polygon=qstab.internal.diagnosticGeometry('inflate',p,radius,o.DiskSides);
region=struct('kind',kind,'polygon',polygon,'basePolygon',p,'controls',controls, ...
    'inflationRadius',radius,'distanceToZero',qstab.internal.diagnosticGeometry('distance',polygon), ...
    'zeroExcludedNumerically',false,'sourceCell',c.cell,'diagnostic',d);
tol=1e-10;if isfield(R,'options')&&isfield(R.options,'Tolerance'),tol=R.options.Tolerance;end
region.zeroExcludedNumerically=region.distanceToZero>tol*max([1;abs(polygon(:))]);
E=struct('schema','qstab.value-enclosure.1','lambda',z,'uncertainty',R.uncertainty, ...
    'degree',R.degree,'regions',{{region}},'scope','oneStoredParameterFrequencyCell', ...
    'parameterCoverage','oneRetainedCellOnly','isCertificate',false,'validated',false, ...
    'assurance','ordinaryDoubleReconstructedGeometry', ...
    'zeroExcludedNumerically',region.zeroExcludedNumerically, ...
    'originalVariable',D.variable,'variable',D.testVariable, ...
    'frequencyBounds',f,'sourceCellIndex',o.CellIndex,'reportedStatus',R.status, ...
    'reportedComplete',D.reportedComplete,'inputBinding',binding,'geometryResidual',residual, ...
    'description',['One saved numerical enclosure, not a whole-family proof. ', ...
      'The displayed disk approximation is outer/circumscribed. ', ...
      'Source report data are read-only and not formally revalidated.']);
end
function P=localPolynomial(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
