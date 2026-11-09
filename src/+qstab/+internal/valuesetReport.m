function S = valuesetReport(P, lambda, varargin)
%VALUESETREPORT Internal implementation returning the full report.

if nargin < 2
    error('qstab:valueset:NotEnoughInputs','Supply P and lambda.');
end
if ~isnumeric(lambda) || ~isvector(lambda) || isempty(lambda) || ...
        any(~isfinite(real(lambda(:)))) || any(~isfinite(imag(lambda(:))))
    error('qstab:valueset:Lambda', ...
        'lambda must be a nonempty finite real or complex numeric vector.');
end
lambda = double(lambda(:).');
if mod(numel(varargin),2) ~= 0
    error('qstab:valueset:NameValue','Options must be name/value pairs.');
end
W = [];
hasWeights = false;
for k = 1:2:numel(varargin)
    name = local_option_name(varargin{k});
    switch name
        case 'weights'
            W = varargin{k+1};
            hasWeights = true;
        otherwise
            error('qstab:valueset:Option','Unknown option "%s".',name);
    end
end

if iscell(P) && ~hasWeights && all(cellfun(@local_scalar_coefficient,P(:)))
    vertices = {P};
elseif iscell(P)
    if isempty(P) || ~isvector(P)
        error('qstab:valueset:Vertices','Vertex polynomials must be a nonempty cell vector.');
    end
    vertices = P(:);
else
    vertices = {P};
end
nv = numel(vertices);
C = cell(nv,1);
degrees = zeros(nv,1);
for j = 1:nv
    C{j} = qstab.polycomponents(vertices{j});
    if ~isnumeric(C{j}) || ~isreal(C{j}) || size(C{j},1) ~= 4 || ...
            isempty(C{j}) || ndims(C{j}) > 2 || any(~isfinite(C{j}(:)))
        error('qstab:valueset:Coefficients', ...
            'Each vertex must yield a finite real 4-by-(grade+1) component array.');
    end
    C{j} = double(C{j});
    degrees(j) = size(C{j},2)-1;
end
grade = max(degrees);
packed = zeros(4*(grade+1),nv);
for j = 1:nv
    cj = zeros(4,grade+1);
    cj(:,1:size(C{j},2)) = C{j};
    packed(:,j) = cj(:);
end

defaultWeights = isempty(W);
if defaultWeights
    W = eye(nv);
end
if ~isnumeric(W) || ~isreal(W) || isempty(W) || ndims(W) > 2 || ...
        any(~isfinite(W(:)))
    error('qstab:valueset:Weights','Weights must be a finite real numeric matrix.');
end
if nv > 1 && isvector(W) && numel(W) == nv
    W = reshape(W,1,nv);
end
if size(W,2) ~= nv
    error('qstab:valueset:WeightSize', ...
        'Weights must have one column per vertex; each row is one sample.');
end
Winput = double(W);
tol = 100*eps*max(1,nv);
sums = sum(Winput,2);
if any(Winput(:) < -tol) || any(abs(sums-1) > tol)
    error('qstab:valueset:Simplex', ...
        'Each weight row must be nonnegative and sum to one within floating-point tolerance.');
end
W = max(Winput,0);
W = bsxfun(@rdivide,W,sum(W,2));

mixed = packed*W.';
ns = size(W,1);
values = complex(zeros(ns,numel(lambda)));
for k = 1:ns
    ck = reshape(mixed(:,k),4,grade+1);
    evaluated = complex(zeros(4,numel(lambda)));
    for ell = grade+1:-1:1
        evaluated = bsxfun(@plus,bsxfun(@times,evaluated,lambda),ck(:,ell));
    end
    values(k,:) = sum(evaluated.^2,1);
end

S = struct();
S.points = values;
S.lambda = lambda;
S.parameters = W;
S.inputWeights = Winput;
S.weightAdjustment = max(abs(W(:)-Winput(:)));
S.vertexCount = nv;
S.sampleCount = ns;
S.vertexGrades = degrees;
S.grade = grade;
S.componentOrder = '1,i,j,k';
S.kind = 'samples';
S.isCertificate = false;
S.finiteMask = isfinite(real(values)) & isfinite(imag(values));
if nv == 1
    S.kind = 'singleton';
    S.description = 'Norm-polynomial values of one fixed scalar polynomial.';
elseif defaultWeights
    S.description = 'Vertex samples only; intermediate convex combinations have not been covered.';
else
    S.description = 'Norm-polynomial values at the supplied shared convex-parameter samples.';
end
end

function name = local_option_name(x)
if isstring(x) && isscalar(x)
    x = char(x);
end
if ~ischar(x) || size(x,1) ~= 1
    error('qstab:valueset:OptionName','Option names must be character rows or scalar strings.');
end
name = lower(x);
end

function tf = local_scalar_coefficient(x)
if isnumeric(x) || islogical(x) || isa(x,'quaternion')
    tf = isscalar(x);
elseif isstruct(x) && isscalar(x) && all(isfield(x,{'a','b','c','d'}))
    tf = isscalar(x.a) && isscalar(x.b) && isscalar(x.c) && isscalar(x.d);
else
    tf = false;
end
end
