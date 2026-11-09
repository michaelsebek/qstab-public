function [D,args,domain]=familyMargin(F,args,operation)
%FAMILYMARGIN Bind the variable, grade, parameter side and parameter domain.
[D,kind,v,g]=qstab.internal.familyDefinition(F);
if ~strcmp(kind,'oneparameter'),error('qstab:family:MarginKind','This margin requires a oneparameter family.');end
expected={'Variable',v,'Side',D.side,'Grade',g};req={'Variable','Side'};
if ~strcmp(v,'s'),req{end+1}='Grade';end
if strcmp(operation,'raymargin'),expected=[expected,{'Type','polynomial'}];end
args=qstab.internal.familyOptions(F,args,expected,req);
extended=false;
if ~strcmp(operation,'raymargin')
 idx=find(cellfun(@(x)((ischar(x)&&isrow(x))||(isstring(x)&&isscalar(x)))&&strcmpi(x,'ExtendParameter'),args(1:2:end)));
 if ~isempty(idx)
  if numel(idx)~=1,error('qstab:family:ExtendParameter','Supply ExtendParameter once.');end
  val=args{2*idx};
  if ~(isscalar(val)&&(islogical(val)||(isnumeric(val)&&isreal(val)&&any(val==[0 1])))),error('qstab:family:ExtendParameter','ExtendParameter must be a logical scalar.');end
  extended=logical(val);args(2*idx-1:2*idx)=[];
 end
 if ~isempty(D.bounds)&&~extended,error('qstab:family:BoundedMargin','This margin searches all quaternion parameters. Supply an unrestricted family or explicitly use ExtendParameter=true.');end
end
domain=struct('scope','unrestrictedModel','side',D.side,'isBounded',~isempty(D.bounds),'bounds',D.bounds,'explicitExtension',extended,'variable',v,'grade',g,'samplingPeriod',F.samplingPeriod);
end
