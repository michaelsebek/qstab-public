function W = witnessDataReport(data,R,varargin)
%WITNESSDATAREPORT Internal implementation returning the full report.
ip=inputParser;ip.FunctionName='qstab.witnessData';
addParameter(ip,'Uncertainty','',@(x)(ischar(x)&&(isrow(x)||isempty(x)))||(isstring(x)&&isscalar(x)));
addParameter(ip,'Tolerance',1e-10,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<.01);
parse(ip,varargin{:});o=ip.Results;
if ~isstruct(R)||~isscalar(R)||~all(isfield(R,{'status','variable','method','witness'}))
    error('qstab:witnessData:Report','Supply a scalar robust-test report.');
end
variable=char(R.variable);
if ~ismember(variable,{'s','z','zinv'}),error('qstab:witnessData:Variable','Unsupported report variable.');end
tv=variable;if strcmp(variable,'zinv'),tv='z';end
W=struct('schema','qstab.witness-data.1','available',~isempty(R.witness), ...
    'reportedStatus',R.status,'checkStatus','noWitness','isAdmissible',false, ...
    'isConfirmedCounterexample',false,'isDegreeLoss',false,'variable',variable, ...
    'testVariable',tv,'components',[],'testComponents',[],'grade',[], ...
    'weights',[],'parameters',[],'membershipResidual',NaN,'parameterResidual',NaN, ...
    'tolerance',o.Tolerance,'nominal',[],'inputRepresentatives',[], ...
    'testRepresentatives',[],'spectrumValid',false,'validated',false,'message','');
if ~W.available,W.message='This report contains no actual member witness.';return;end
kind=char(o.Uncertainty);
if isempty(kind)&&isfield(R,'uncertainty'),kind=R.uncertainty;end
F=qstab.familycomponents(data,'Uncertainty',kind);w=R.witness;
if ~isstruct(w)||~isscalar(w),error('qstab:witnessData:Witness','Witness must be a scalar struct.');end
W.grade=F.degree;C=[];
if isfield(w,'originalComponents')
    C=w.originalComponents;
elseif isfield(w,'coefficients')
    C=qstab.polycomponents(w.coefficients);
elseif isfield(w,'components')
    C=w.components;
    if strcmp(variable,'zinv')&&~strcmp(R.method,'subdivisionLyapunovPerturbation'),C=C(:,end:-1:1);end
end
W.checkStatus='invalidWitness';
if ~isnumeric(C)||~isreal(C)||~isequal(size(C),[4 F.degree+1])||any(~isfinite(C(:)))
    W.message='Witness is absent, nonfinite or has the wrong original declared grade.';return
end
W.components=C;W.testComponents=C;
if strcmp(variable,'zinv'),W.testComponents=C(:,end:-1:1);end
factor=100*o.Tolerance;scale=max(1,max(abs(C(:))));
switch F.kind
    case 'interval'
        W.membershipResidual=max([0;F.intervalLower(:)-C(:);C(:)-F.intervalUpper(:)]);
        W.parameterResidual=W.membershipResidual;
        scale=max([scale;abs(F.intervalLower(:));abs(F.intervalUpper(:))]);
        W.isAdmissible=W.membershipResidual<=factor*scale;
    case 'polytope'
        t=localFirst(w,{'originalWeights','vertexWeights'});t=t(:);W.weights=t;
        if ~isnumeric(t)||~isreal(t)||numel(t)~=size(F.vertices,3)||any(~isfinite(t))
            W.message='Original convex weights are missing or invalid.';return
        end
        D=zeros(size(C));for j=1:numel(t),D=D+t(j)*F.vertices(:,:,j);end
        W.membershipResidual=max(abs(D(:)-C(:)));
        W.parameterResidual=max([0;-t;abs(sum(t)-1)]);
        scale=max([scale;abs(D(:))]);
        W.isAdmissible=W.parameterResidual<=factor&&W.membershipResidual<=factor*scale;
    otherwise
        t=localFirst(w,{'parameter','theta'});t=t(:);W.parameters=t;
        if ~isnumeric(t)||~isreal(t)||numel(t)~=size(F.bounds,1)||any(~isfinite(t))
            W.message='Original shared parameter vector is missing or invalid.';return
        end
        D=F.base;for j=1:numel(t),D=D+t(j)*F.generators(:,:,j);end
        W.membershipResidual=max(abs(D(:)-C(:)));
        ps=max(1,max(abs(F.bounds),[],2));
        W.parameterResidual=max([0;(F.bounds(:,1)-t)./ps;(t-F.bounds(:,2))./ps]);
        scale=max([scale;abs(D(:))]);
        W.isAdmissible=W.parameterResidual<=factor&&W.membershipResidual<=factor*scale;
end
if ~W.isAdmissible,W.message='The supplied coefficients do not reconstruct as an admissible original-family member.';return;end
lead=F.degree+1;if strcmp(variable,'zinv'),lead=1;end
W.isDegreeLoss=all(C(:,lead)==0);
if strcmp(R.status,'assumptionFailed')
    if W.isDegreeLoss
        W.checkStatus='confirmedDegreeLoss';W.message='An admissible exactly zero required coefficient; not a finite-root instability claim.';
    else
        W.checkStatus='unconfirmedAssumption';W.message='Membership passed, but the required coefficient is not exactly zero.';
    end
    return
end
W.nominal=qstab.internal.isstableReport(localPolynomial(C),'Variable',variable,'Method','spectral','Tolerance',o.Tolerance);
if isfield(W.nominal,'diagnostics')&&isfield(W.nominal.diagnostics,'representatives')
    W.testRepresentatives=W.nominal.diagnostics.representatives;
    W.inputRepresentatives=W.testRepresentatives;
    if isfield(W.nominal.diagnostics,'spectrum')
        W.spectrumValid=isequal(W.nominal.diagnostics.spectrum.isValid,true);
    end
    if strcmp(variable,'zinv')
        z=W.testRepresentatives;d=complex(zeros(size(z)));zero=(z==0);infinite=isinf(abs(z));
        good=~zero&~infinite;d(good)=1./z(good);
        d(good)=complex(real(d(good)),abs(imag(d(good))));d(zero)=Inf;d(infinite)=0;
        W.inputRepresentatives=d;
    end
end
W.isConfirmedCounterexample=strcmp(R.status,'notStable')&&strcmp(W.nominal.status,'notStable')&&W.spectrumValid;
if W.isConfirmedCounterexample
    W.checkStatus='confirmedNotStable';W.message='Numerically admissible original member; independent nominal spectral check agrees.';
else
    W.checkStatus='unresolvedSpectralCheck';W.message=['Membership passed; independent nominal status: ',W.nominal.status,'.'];
end
end
function v=localFirst(s,names)
v=[];for j=1:numel(names),if isfield(s,names{j})&&~isempty(s.(names{j})),v=s.(names{j});return;end;end
end
function P=localPolynomial(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
