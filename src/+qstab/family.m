function F = family(kind,varargin)
%FAMILY Define a scalar quaternion polynomial family using ordinary MATLAB data.
%   F = qstab.family('interval',L,U) uses independent bounds on each real,
%   i, j and k part of each coefficient. L must be <= U componentwise.
%   F = qstab.family('polytope',{P1,...,Pm}) uses one real nonnegative weight
%   per generator, summing to one. The weights are common to all powers.
%   F = qstab.family('affine',P0,{G1,...,Gm},B) denotes P0+sum(t_j*Gj);
%   each real t_j belongs to [B(j,1),B(j,2)], shared across coefficients.
%   F = qstab.family('oneparameter',P0,P1,'Side','left','Bounds',B)
%   denotes P0+t*P1, with ONE quaternion t. 'right' denotes P0+P1*t.
%   B is a real ordered 4-by-2 box on real(t),i,j,k coordinates. Omit B
%   or use [] for all of H: valid for a ray margin; a bounded-domain test
%   or plot requires an explicit finite domain.
%
%   All input polynomials have ASCENDING coefficients; use MATLAB quaternion,
%   numeric real/complex vectors or component structs a,b,c,d. 'Variable'
%   is 's' (default), 'z', or 'zinv'. 'Grade' pads with high-order zero
%   coefficients; default = maximum supplied grade. 'SamplingPeriod' is 0
%   for 's' and positive for discrete variables (default 1). It is metadata;
%   no transformation or frequency scaling is applied.
%   F is a scalar STRUCT, not a toolbox class. Fields include kind, variable,
%   grade, samplingPeriod and lower/upper, vertices, center/generators/bounds,
%   or p0/p1/bounds/side according to kind. No stability test is performed.
%   F = qstab.family('polytope',{[1 1],[2 1]}); tf=qstab.robuststable(F);
%   F = qstab.family('oneparameter',[2 1],1,'Bounds',[-1 .5;zeros(3,2)]);
%   [ab,I]=qstab.raymargin(F,-2); I.parameterDomain.admissibleAlpha
%   See also qstab.randfamily, qstab.familycomponents, qstab.robuststable.
try
if nargin==0 || ~textScalar(kind), error('qstab:family:Kind','Supply a family kind.'); end
kind=lower(char(kind));
F=struct('schema','qstab.family-definition.1','kind',kind,'variable','s','grade',[],'samplingPeriod',[]);
switch kind
 case {'interval','oneparameter'}, n=2;
 case 'polytope', n=1;
 case 'affine', n=3;
 otherwise, error('qstab:family:Kind','Use interval, polytope, affine or oneparameter.');
end
if numel(varargin)<n || rem(numel(varargin)-n,2), error('qstab:family:Arguments','Supply defining data followed by complete name/value pairs.'); end
switch kind
 case 'interval', F.lower=varargin{1};F.upper=varargin{2};polys=varargin(1:2);
 case 'polytope'
    V=varargin{1};
    if ~iscell(V)||~isvector(V)||isempty(V),error('qstab:family:Vertices','Supply a nonempty cell vector of polynomial generators.');end
    F.vertices=V;polys=V(:).';
 case 'affine'
    G=varargin{2};if ~iscell(G)||(~isempty(G)&&~isvector(G)),error('qstab:family:Generators','Supply a cell vector of polynomial directions.');end
    F.center=varargin{1};F.generators=G;F.bounds=varargin{3};polys=[varargin(1) G(:).'];
 case 'oneparameter',F.p0=varargin{1};F.p1=varargin{2};F.bounds=[];F.side='left';polys=varargin(1:2);
end
seen={};for i=n+1:2:numel(varargin)
 if ~textScalar(varargin{i}),error('qstab:family:Option','A name/value option requires a text name.');end
 key=lower(char(varargin{i}));v=varargin{i+1};
 if any(strcmp(key,seen)),error('qstab:family:Option','Duplicate option: %s.',key);end
 seen{end+1}=key; %#ok<AGROW>
 switch key
  case 'variable',F.variable=v;
  case 'grade',F.grade=v;
  case 'samplingperiod',F.samplingPeriod=v;
  case {'side','bounds'}
   if ~strcmp(kind,'oneparameter'),error('qstab:family:Option','Use %s as an option only for oneparameter.',key);end
   F.(key)=v;
  otherwise,error('qstab:family:Option','Unknown option: %s.',key);
 end
end
if isempty(F.grade)
 g=0;for k=1:numel(polys),c=qstab.polycomponents(polys{k});g=max(g,size(c,2)-1);end
 F.grade=g;
end
if ~textScalar(F.variable),error('qstab:family:Variable','Variable must be s, z, or zinv.');end
F.variable=lower(char(F.variable));
if isempty(F.samplingPeriod),F.samplingPeriod=double(~strcmp(F.variable,'s'));end
qstab.internal.familyDefinition(F); % validate definition without creating a dense box expansion
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
function b=textScalar(a)
b=(ischar(a)&&isrow(a))||(isstring(a)&&isscalar(a)&&~ismissing(a));
end
