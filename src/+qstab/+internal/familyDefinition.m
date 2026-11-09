function [D,kind,v,g]=familyDefinition(F)
%FAMILYDEFINITION Validate a public definition and expose ordinary input data.
if ~qstab.internal.isFamilyDefinition(F)||~all(isfield(F,{'kind','variable','grade','samplingPeriod'}))
 error('qstab:family:Definition','Expected the complete scalar struct produced by qstab.family.');
end
kind=F.kind;v=F.variable;g=F.grade;h=F.samplingPeriod;
if ~ischar(kind)||~isrow(kind),error('qstab:family:Kind','Invalid family kind.');end
if ~ischar(v)||~isrow(v)||~any(strcmp(v,{'s','z','zinv'})),error('qstab:family:Variable','Variable must be s, z or zinv.');end
if ~isnumeric(g)||~isreal(g)||~isscalar(g)||~isfinite(g)||g<0||g~=fix(g),error('qstab:family:Grade','Grade must be a nonnegative integer.');end
if ~(isnumeric(h)&&isreal(h)&&isscalar(h)&&isfinite(h)&&((strcmp(v,'s')&&h==0)||(~strcmp(v,'s')&&h>0))),error('qstab:family:SamplingPeriod','Use zero for continuous time, a positive period for discrete time.');end
switch kind
 case 'interval'
  require(F,{'lower','upper'});[L,l]=poly(F.lower,g);[U,u]=poly(F.upper,g);
  if any(l(:)>u(:)),error('qstab:family:IntervalOrder','Each real coefficient component of lower must be <= upper.');end
  D={L,U};
 case 'polytope'
  require(F,{'vertices'});D=F.vertices;
  if ~iscell(D)||~isvector(D)||isempty(D),error('qstab:family:Vertices','Expected a nonempty cell vector of polynomials.');end
  for i=1:numel(D),D{i}=poly(D{i},g);end
 case 'affine'
  require(F,{'center','generators','bounds'});G=F.generators;
  if ~iscell(G)||(~isempty(G)&&~isvector(G)),error('qstab:family:Generators','Generators must be a cell vector.');end
  for i=1:numel(G),G{i}=poly(G{i},g);end
  D=struct('center',poly(F.center,g),'generators',{G},'bounds',bounds(F.bounds,numel(G)));
 case 'oneparameter'
  require(F,{'p0','p1','bounds','side'});B=F.bounds;
  if ~(isnumeric(B)&&isreal(B)&&isequal(size(B),[0 0])),B=bounds(B,4);end
  if ~((ischar(F.side)&&isrow(F.side))||(isstring(F.side)&&isscalar(F.side)&&~ismissing(F.side)))||~any(strcmpi(F.side,{'left','right'})),error('qstab:family:Side','Side must be left or right.');end
  D=struct('p0',poly(F.p0,g),'p1',poly(F.p1,g),'bounds',B,'side',lower(char(F.side)));
 otherwise,error('qstab:family:Kind','Invalid family kind.');
end
end
function [P,C]=poly(P,g)
C=qstab.polycomponents(P,'Grade',g);
if any(~isfinite(C(:))),error('qstab:family:Coefficients','Coefficients must be finite.');end
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
function B=bounds(B,n)
if ~isnumeric(B)||~isreal(B)||~isequal(size(B),[n 2])||any(~isfinite(B(:)))||any(B(:,1)>B(:,2)),error('qstab:family:Bounds','Supply %d rows of finite ordered real endpoints.',n);end
B=double(B);
end
function require(S,fs)
if ~all(isfield(S,fs)), error('qstab:family:MissingField','The definition is missing a required field.'); end
end
