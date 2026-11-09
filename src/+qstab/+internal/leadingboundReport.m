function R = leadingboundReport(data,varargin)
%LEADINGBOUNDREPORT Internal implementation returning the full report.
ip=inputParser;ip.FunctionName='qstab.leadingbound';
addParameter(ip,'Uncertainty','',@(x)(ischar(x)&&(isrow(x)||isempty(x)))||(isstring(x)&&isscalar(x)));
addParameter(ip,'Variable','s',@(x)(ischar(x)&&isrow(x))||(isstring(x)&&isscalar(x)));
addParameter(ip,'Tolerance',1e-10,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<.01);
addParameter(ip,'MaxIterations',256,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x==floor(x));
addParameter(ip,'MaxTime',30,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&x>=0&&~isnan(x));
parse(ip,varargin{:});o=ip.Results;
v=lower(strtrim(char(o.Variable)));
if ismember(v,{'continuous','hurwitz'}),v='s';end
if ismember(v,{'discrete','schur'}),v='z';end
if ismember(v,{'delay','z^-1','z^{-1}','d'}),v='zinv';end
if ~ismember(v,{'s','z','zinv'})
    error('qstab:leadingbound:Variable','Variable must be s, z, or zinv.');
end
F=qstab.familycomponents(data,'Uncertainty',o.Uncertainty);
R=qstab.internal.leadbound(F,v,o.Tolerance,o.MaxIterations,o.MaxTime);
R.uncertainty=F.kind;R.degree=F.degree;R.options=o;
end
