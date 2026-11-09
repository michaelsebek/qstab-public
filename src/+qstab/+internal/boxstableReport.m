function R = boxstableReport(P,widths,r,varargin)
%BOXSTABLEREPORT Internal implementation returning the full report.
ip=inputParser;
addParameter(ip,'Variable','s',@localText);addParameter(ip,'Grade',[]);
addParameter(ip,'Method','zeroexclusion',@localText);
addParameter(ip,'MaxCells',2048,@(x)localInteger(x)&&x>=1);
addParameter(ip,'MaxDepth',30,@localInteger);
addParameter(ip,'MaxTime',20,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&~isnan(x)&&x>=0);
addParameter(ip,'Tolerance',1e-10,@(x)localPositive(x)&&x<1);
addParameter(ip,'StoreCertificates',false,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});o=ip.Results;
numericNames={'AbsTol','RelTol','Tolerance','InitialRadius','MaxIterations','MaxCells','MaxDepth','MaxTime','MaxRadius','MaxFacets'};
for jj=1:numel(numericNames)
    if isfield(o,numericNames{jj}),o.(numericNames{jj})=double(o.(numericNames{jj}));end
end
v=lower(char(o.Variable));method=lower(char(o.Method));
if ~any(strcmp(v,{'s','z','zinv'})),error('qstab:boxstable:Variable','Variable must be s, z or zinv.');end
if ~any(strcmp(method,{'zeroexclusion','lyapunov'})),error('qstab:boxstable:Method','Method must be zeroexclusion or lyapunov.');end
if ~(isnumeric(r)&&isreal(r)&&isscalar(r)&&isfinite(r)&&r>=0)
    error('qstab:boxstable:Radius','Radius must be finite, real and nonnegative.');end
r=double(r);
M=qstab.internal.boxmodel(P,widths,o.Grade);
if M.grade>64,error('qstab:boxstable:DegreeLimit','Current dense whole-box tests support grade at most 64.');end
N=qstab.internal.isstableReport(localParts(M.C),'Variable',v,'Grade',M.grade,'Tolerance',o.Tolerance);
R=struct('schema','qstab.boxstable.1','status','undetermined','isStable',[], ...
    'variable',v,'grade',M.grade,'radius',double(r),'weights',M.weights,'nominal',N, ...
    'method',method,'witness',[],'familyResult',[],'diagnostics',struct(), ...
    'message','','validated',false);
if strcmp(N.status,'notStable')
    R.status='notStable';R.isStable=false;R.witness=struct('kind','nominal','valid',true, ...
        'weightedNorm',0,'coefficientComponents',M.C,'polynomial',localParts(M.C));
    R.message='The nominal member belongs to every scaled box and is unstable.';return
elseif ~strcmp(N.status,'stable')
    R.status=N.status;R.message='The nominal member is not verified strictly stable.';return
end
if r==0||~M.hasPerturbation
    R.status='stable';R.isStable=true;R.message='The requested family contains only the verified nominal member.';return
end
[lead,leadWitness]=qstab.internal.boxdegree(M,v);
if isfinite(lead)&&r>=lead
    R.status='assumptionFailed';R.witness=leadWitness;R.message='The box loses its required leading coefficient.';return
end
T=qstab.internal.boxtrial(M,r,v,method,o.MaxCells,o.MaxDepth,o.MaxTime,o.Tolerance,o.StoreCertificates);
R.familyResult=T.result;R.diagnostics=T;R.message=T.message;
if strcmp(T.status,'stable')
    R.status='stable';R.isStable=true;
elseif ~isempty(T.witness)&&T.witness.valid&&T.witness.insideRequestedBox&&T.witness.weightedNorm<=r
    R.status=T.status;R.witness=T.witness;
    if strcmp(T.status,'notStable'),R.isStable=false;end
end
end
function P=localParts(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
function y=localText(x),y=ischar(x)||(isstring(x)&&isscalar(x));end
function y=localPositive(x),y=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0;end
function y=localInteger(x),y=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x==fix(x);end
