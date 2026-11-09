function F = randfamily(varargin)
%RANDFAMILY Draw a repeatable, well-defined scalar quaternion polynomial family.
%   F = qstab.randfamily                         random interval, degree 1, 2 or 3
%   F = qstab.randfamily(D)                      interval of degree D
%   F = qstab.randfamily(D,KIND)                 choose one of the four family kinds
%   F = qstab.randfamily(D,KIND,N)               N polynomials or affine directions
%   F = qstab.randfamily(KIND)                   choose the kind; draw the degree
%   F = qstab.randfamily(KIND,D,...)             equivalent kind-first spelling
%   F = qstab.randfamily(...,'int',...,'Seed',42) small integer coefficients, repeatable
%
%   D is the degree in the forward variable s or z. [] chooses uniformly
%   from 1:3. In delay notation the coefficients are reversed at grade D.
%   D=0 generates nonzero quaternion constants, except for 'oneparameter',
%   which requires a nonzero perturbation of strictly lower forward degree.
%
%   KIND           Family and meaning of N
%   'interval'     Default. Independent bounds on the real,i,j,k parts of
%                  every free coefficient. No N.
%   'oneparameter' One quaternion t shared by all powers. P0+t*P1 (left)
%                  or P0+P1*t (right). No N. The nominal is built to be stable.
%   'polytope'     Real convex hull of N defining polynomials (default 3).
%                  N counts generators, not necessarily extreme points.
%   'affine'       P0+sum(theta_j*G_j) with N shared real parameters
%                  (default 2). One interval per parameter.
%
%   Optional inputs
%   's','z','d' or 'Variable',V  Time variable. Aliases: p for s, q for z,
%        zi/z^-1/zinv for delay. Default: s.
%   'Width',W    Positive coefficient scale (default 1). Scales the interval
%        half-widths, the displacements from the first convex generator, or
%        the perturbation directions. It is not a stability radius.
%   'int'        Integer defining coefficients; W must be a positive integer.
%        This does not discretize the parameters, weights or family members.
%   'Side',S     'left' (default) or 'right'; only for 'oneparameter'.
%   'Bounds',B   'oneparameter': finite ordered 4-by-2 bounds on real,i,j,k
%        of t, default [-1,1] for each. [] means all of H; use for margins.
%        'affine': finite ordered N-by-2 bounds; each defaults to [-1,1].
%        Equal endpoints fix that coordinate. Bounds are independent of W.
%   'SamplingPeriod',H   0 in continuous time; positive in discrete time
%        (discrete default 1, in the same units as your model).
%   'Seed',K     Integer 0..2^32-1. Repeats the family without advancing
%        the caller's random stream. Without K the active stream advances.
%   'Degree',D, 'Kind',KIND (or 'Type',KIND) and 'Count',N are also accepted.
%        Use each argument once. N is a positive integer for polytope/affine.
%
%   F is an ordinary definition struct accepted by qstab.robuststable,
%   qstab.raymargin (oneparameter), and qstab.plotFamilyValues / plotFamilyZeros.
%   Its defining polynomials are component structs with row vectors a,b,c,d.
%   For D>=1 every member is monic in the forward variable. Delay coefficients
%   are reversed at the full grade; the current-sample coefficient is one.
%   The nominal is stable by construction in exact arithmetic. No stability
%   test of a whole random family is implied.
%
%   F=qstab.randfamily(3,'interval','Seed',17);
%   [~,V]=qstab.plotFamilyValues(F,1,9,'Plot',false);
%   F=qstab.randfamily(2,'polytope',4,'z','Seed',2);
%   [tf,I]=qstab.robuststable(F,'MaxCells',16);
%   F=qstab.randfamily(2,'oneparameter','Bounds',[],'Seed',17);
%   ab=qstab.raymargin(F,1);
%   See also qstab.family, qstab.familycomponents, qstab.plotFamilyZeros.
try
nargoutchk(0,1);
o = parseInputs(varargin);
% With a seed use a separate stream, including for the optional degree.
% Never reset or replace the caller's MATLAB RandStream object.
o.uniform = @rand; o.drawint = @randi;
if ~isempty(o.seed)
    if exist('RandStream','class') == 8
        stream = RandStream('mt19937ar','Seed',o.seed);
        o.uniform = @(varargin) rand(stream,varargin{:});
        o.drawint = @(varargin) randi(stream,varargin{:});
    else
        % Hosts without RandStream (e.g. for exercising numeric code in Octave).
        saved = rng; guard = onCleanup(@()rng(saved)); %#ok<NASGU>
        rng(o.seed,'twister');
    end
end
if isempty(o.degree), o.degree = o.drawint(3); end
D=qstab.internal.randomFamilyCoefficients(o);
meta={'Variable',o.variable,'Grade',o.degree,'SamplingPeriod',o.samplingPeriod};
switch o.kind
 case 'interval',F=qstab.family('interval',asPol(D.lower),asPol(D.upper),meta{:});
 case 'polytope',V=cellfun(@asPol,D.vertices,'UniformOutput',false);F=qstab.family('polytope',V,meta{:});
 case 'affine',G=cellfun(@asPol,D.generators,'UniformOutput',false);F=qstab.family('affine',asPol(D.center),G,D.bounds,meta{:});
 case 'oneparameter',F=qstab.family('oneparameter',asPol(D.p0),asPol(D.p1),'Bounds',D.bounds,'Side',D.side,meta{:});
end
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end

function P=asPol(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end

function o = parseInputs(args)
o = struct('degree',[],'kind','interval','count',[],'variable','','symbol','', ...
    'side','left','bounds',[],'width',1,'integer',false,'samplingPeriod',[],'seed',[]);
seen = struct;
k = 1;
while k <= numel(args)
    a = args{k}; k = k+1;
    if isnumeric(a)
        if ~isfield(seen,'degree'), key = 'degree';
        elseif ~isfield(seen,'count'), key = 'count';
        else, error('qstab:randfamily:Positional','Only degree and (for polytope/affine) count are positional numbers.'); end
        value = a;
    elseif isText(a)
        t = lower(strtrim(char(a)));
        switch t
            case {'interval','polytope','affine','oneparameter'}
                key='kind'; value=t;
            case {'int','integer'}
                key='integer'; value=true;
            case {'s','p','z','q','d','zi','z^-1','zinv'}
                key='symbol'; value=t;
            case {'degree','kind','type','count','variable','width','side','bounds','samplingperiod','seed'}
                key=t;
                if strcmp(key,'type'),key='kind';end
                if strcmp(key,'variable'),key='symbol';end
                if strcmp(key,'samplingperiod'),key='samplingPeriod';end
                if k>numel(args),error('qstab:randfamily:MissingValue','Missing value after "%s".',t);end
                value=args{k}; k=k+1;
            otherwise
                error('qstab:randfamily:Option','Unknown option "%s". See help qstab.randfamily.',t);
        end
    else
        error('qstab:randfamily:Argument','Expected a numeric degree/count, a family kind, or an option.');
    end
    if isfield(seen,key),error('qstab:randfamily:Duplicate','"%s" was supplied twice.',key);end
    seen.(key)=true; o.(key)=value;
end
o.kind = choice(o.kind,{'interval','oneparameter','polytope','affine'},'Kind');
if ~isempty(o.degree) || ~isnumeric(o.degree)
    o.degree=whole(o.degree,0,flintmax-1,'Degree');
end
if strcmp(o.kind,'oneparameter') && isequal(o.degree,0)
    error('qstab:randfamily:Degree','Use degree >= 1 for a nonzero lower-degree oneparameter perturbation.');
end
many=any(strcmp(o.kind,{'polytope','affine'}));
if ~many && (~isnumeric(o.count) || ~isempty(o.count))
    error('qstab:randfamily:Count','Count applies only to polytope (generators) or affine (directions), not %s.',o.kind);
elseif many
    if isnumeric(o.count)&&isempty(o.count), o.count=2+strcmp(o.kind,'polytope'); end
    o.count=whole(o.count,1,flintmax-1,'Count');
end
if isfield(seen,'side') && ~strcmp(o.kind,'oneparameter')
    error('qstab:randfamily:Side','Side is meaningful only for oneparameter.');
end
o.side=choice(o.side,{'left','right'},'Side');
if isfield(seen,'bounds') && ~any(strcmp(o.kind,{'oneparameter','affine'}))
    error('qstab:randfamily:Bounds','Bounds applies only to affine and oneparameter.');
end
if any(strcmp(o.kind,{'oneparameter','affine'}))
    n=o.count; if strcmp(o.kind,'oneparameter'), n=4; end
    if ~isfield(seen,'bounds')
        o.bounds = repmat([-1 1],n,1);
    elseif ~(strcmp(o.kind,'oneparameter') && isnumeric(o.bounds) && isreal(o.bounds) && isequal(size(o.bounds),[0 0]))
        b=o.bounds;
        if ~(isnumeric(b)&&isreal(b)&&isequal(size(b),[n 2])&&all(isfinite(b(:)))&&all(b(:,1)<=b(:,2)))
            error('qstab:randfamily:Bounds','Expected %d rows [lower upper] of finite real ordered bounds.',n);
        end
        o.bounds=full(double(b));
    end
end
if ~(scalarNumber(o.width) && o.width>0)
    error('qstab:randfamily:Width','Width must be a positive finite real scalar.');
end
o.width=double(o.width);
if o.integer && (fix(o.width)~=o.width || o.width>flintmax)
    error('qstab:randfamily:Width','With int, Width must be a positive exactly representable integer <= flintmax.');
end
symbols={'s','p','z','q','d','zi','z^-1','zinv'};
if ~isfield(seen,'symbol'), o.symbol='s';else,o.symbol=choice(o.symbol,symbols,'Variable');end
if any(strcmp(o.symbol,{'s','p'})), o.variable='s';
elseif any(strcmp(o.symbol,{'z','q'})), o.variable='z';
else,o.variable='zinv';end
if ~isfield(seen,'samplingPeriod')
    o.samplingPeriod=double(~strcmp(o.variable,'s'));
elseif ~(scalarNumber(o.samplingPeriod) && ((strcmp(o.variable,'s') && o.samplingPeriod==0) || ...
        (~strcmp(o.variable,'s') && o.samplingPeriod>0)))
    error('qstab:randfamily:SamplingPeriod','Use 0 for continuous time or a positive sampling period for discrete time.');
end
o.samplingPeriod=double(o.samplingPeriod);
if isfield(seen,'seed'),o.seed=whole(o.seed,0,2^32-1,'Seed');end
end

function b=isText(x)
b=(ischar(x)&&isrow(x)) || (isstring(x)&&isscalar(x)&&~ismissing(x));
end

function t=choice(x,allowed,label)
if isText(x)
    t=lower(strtrim(char(x)));
    if any(strcmp(t,allowed)),return;end
end
error(['qstab:randfamily:' label],'%s must be one of: %s.',label,strjoin(allowed,', '));
end

function b=scalarNumber(x)
b=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x);
end

function x=whole(x,lo,hi,label)
if ~scalarNumber(x) || x<lo || x>hi || x~=fix(x)
    error(['qstab:randfamily:' label],'%s must be an integer in [%g, %g].',label,lo,hi);
end
x=full(double(x));
end
