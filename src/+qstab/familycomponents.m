function F = familycomponents(data, varargin)
%FAMILYCOMPONENTS Real component description of a scalar polynomial family.
%   F = qstab.familycomponents(definition) accepts a finite qstab.family struct.
%   F = qstab.familycomponents(DATA,'Uncertainty',KIND) accepts:
%     'polytope'    DATA = {P1,...,Pm}; the real convex hull of the vertices.
%     'interval'    DATA = {L,U}; componentwise bounds, including all four
%                   real components of every quaternion coefficient.
%     'affine'      DATA.center, DATA.generators (a cell array), and
%                   DATA.bounds (number-of-generators by 2, real).
%     'oneparameter' DATA.p0, DATA.p1, DATA.side ('left' or 'right'), and
%                   DATA.bounds (4 by 2, for real,i,j,k components of t).
%   An omitted KIND is inferred from a struct, or is 'polytope' for a cell.
%
%   All polynomials use ASCENDING coefficients. Each polynomial is accepted
%   by qstab.polycomponents; the preferred public type is MATLAB quaternion.
%   Native structs with real fields a,b,c,d are also accepted. No custom
%   class, quaternion constructor, or external quaternion toolbox is used.
%
%   F is an ordinary struct containing only real double arrays and metadata.
%   F.vertices has size 4-by-(degree+1)-by-m for a polytope. Other kinds use
%       P(theta) = F.base + sum_j theta(j)*F.generators(:,:,j).
%   F.bounds contains the actual real parameters used in this expression.
%   For intervals, these parameters are the varying coefficient components
%   themselves, identified by F.parameterIndex = [component,coefficient].
%   The coefficient index in parameterIndex starts at zero.
%
%   Different vertex/generator lengths are padded with zero HIGH-degree
%   coefficients. The resulting declared degree is retained. This function
%   does not certify a nonzero leading coefficient over the family.
%
%   The oneparameter form is a general affine family, also when p1 changes
%   the leading coefficient. The additional deg(p1)<deg(p0) assumption of a
%   particular stability-radius theorem is not imposed by this converter.
%
%   Examples (ascending coefficients):
%   F = qstab.familycomponents({[1 1],[2 1]},'Uncertainty','polytope');
%   disp(F.kind); disp(F.degree);
%
%   Interpretation:
%   This is an inspection/conversion representation with real arrays, not a new polynomial object. Preserve common weights and parameters. The function does not itself prove robust stability.
%
%   Reference: docs/reference/familycomponents.html.
%   See also qstab.about.

try
    ip = inputParser;
    ip.FunctionName = 'qstab.familycomponents';
    addParameter(ip, 'Uncertainty', '', @localIsText);
    parse(ip, varargin{:});
    kind = lower(strtrim(char(ip.Results.Uncertainty)));
    tagged = qstab.internal.isFamilyDefinition(data);
    if tagged
        taggedData=data;
        [data,k,v,g]=qstab.internal.familyDefinition(taggedData);
        h=taggedData.samplingPeriod;
        if ~isempty(kind) && ~strcmp(k,kind),error('qstab:family:OptionConflict','Uncertainty conflicts with the family kind.');end
        if strcmp(k,'oneparameter') && isempty(data.bounds)
            error('qstab:family:Unbounded','Supply finite Bounds on all four coordinates for this operation; raymargin accepts an unrestricted family.');
        end
        kind=k;
    end
    if isempty(kind)
        if iscell(data)
            kind = 'polytope';
        elseif isstruct(data) && isscalar(data) && ...
                isfield(data, 'p0') && isfield(data, 'p1')
            kind = 'oneparameter';
        elseif isstruct(data) && isscalar(data) && ...
                isfield(data, 'center') && isfield(data, 'generators')
            kind = 'affine';
        else
            error('qstab:familycomponents:UncertaintyRequired', ...
                'Specify Uncertainty as polytope, interval, affine, or oneparameter.');
        end
    end
    if strcmp(kind, 'one-parameter') || strcmp(kind, 'oneparam')
        kind = 'oneparameter';
    end

    F = struct('schema', 'qstab.scalar-family.1', 'kind', kind, ...
        'degree', [], 'base', [], 'generators', [], 'bounds', zeros(0,2), ...
        'vertices', [], 'parameterNames', {{}}, 'parameterIndex', zeros(0,2), ...
        'side', '', 'intervalLower', [], 'intervalUpper', []);

    if tagged,F.variable=v;F.grade=g;F.samplingPeriod=h;end
    switch kind
        case 'polytope'
            if ~iscell(data) || isempty(data) || ~isvector(data)
                error('qstab:familycomponents:InvalidVertices', ...
                    'A polytope must be a nonempty cell vector of scalar polynomials.');
            end
            m = numel(data);
            parts = cell(1,m);
            count = 1;
            for j = 1:m
                parts{j} = localPolynomial(data{j});
                count = max(count, size(parts{j},2));
            end
            F.vertices = zeros(4,count,m);
            for j = 1:m
                F.vertices(:,1:size(parts{j},2),j) = parts{j};
            end
            F.degree = count-1;
            F.originalVertexCount = m;

        case 'interval'
            if ~iscell(data) || numel(data) ~= 2
                error('qstab:familycomponents:InvalidInterval', ...
                    'An interval family must be the two-element cell {L,U}.');
            end
            L = localPolynomial(data{1});
            U = localPolynomial(data{2});
            if ~isequal(size(L),size(U))
                error('qstab:familycomponents:IntervalSize', ...
                    'Lower and upper coefficient arrays must have the same length.');
            end
            if any(L(:)>U(:))
                error('qstab:familycomponents:IntervalOrder', ...
                    'Every real coefficient component of L must be <= that of U.');
            end
            F.degree = size(L,2)-1;
            F.intervalLower = L;
            F.intervalUpper = U;
            F.base = L;
            varying = find(L(:)<U(:));
            count = numel(varying);
            F.generators = zeros(4,size(L,2),count);
            F.bounds = zeros(count,2);
            F.parameterIndex = zeros(count,2);
            F.parameterNames = cell(count,1);
            componentNames = {'real','i','j','k'};
            for j = 1:count
                [component,coefficient] = ind2sub(size(L),varying(j));
                F.base(component,coefficient) = 0;
                F.generators(component,coefficient,j) = 1;
                F.bounds(j,:) = [L(component,coefficient),U(component,coefficient)];
                F.parameterIndex(j,:) = [component,coefficient-1];
                F.parameterNames{j} = sprintf('%s(a_%d)', ...
                    componentNames{component},coefficient-1);
            end

        case 'affine'
            localRequiredStruct(data, {'center','generators','bounds'}, 'affine');
            if ~iscell(data.generators) || ...
                    (~isempty(data.generators) && ~isvector(data.generators))
                error('qstab:familycomponents:InvalidGenerators', ...
                    'generators must be a cell vector, possibly empty.');
            end
            base = localPolynomial(data.center);
            m = numel(data.generators);
            B = localBounds(data.bounds,m);
            parts = cell(1,m);
            count = size(base,2);
            for j = 1:m
                parts{j} = localPolynomial(data.generators{j});
                count = max(count,size(parts{j},2));
            end
            F.base = zeros(4,count);
            F.base(:,1:size(base,2)) = base;
            F.generators = zeros(4,count,m);
            for j = 1:m
                F.generators(:,1:size(parts{j},2),j) = parts{j};
            end
            F.bounds = B;
            F.degree = count-1;
            F.parameterNames = localParameterNames(m,'theta');

        case 'oneparameter'
            localRequiredStruct(data, {'p0','p1','side','bounds'}, 'oneparameter');
            if ~localIsText(data.side)
                error('qstab:familycomponents:InvalidSide', ...
                    'side must be left or right.');
            end
            side = lower(strtrim(char(data.side)));
            if strcmp(side,'l'), side = 'left'; end
            if strcmp(side,'r'), side = 'right'; end
            if ~ismember(side,{'left','right'})
                error('qstab:familycomponents:InvalidSide', ...
                    'side must be left or right.');
            end
            p0 = localPolynomial(data.p0);
            p1 = localPolynomial(data.p1);
            count = max(size(p0,2),size(p1,2));
            F.base = zeros(4,count);
            F.base(:,1:size(p0,2)) = p0;
            P1 = zeros(4,count);
            P1(:,1:size(p1,2)) = p1;
            F.generators = zeros(4,count,4);
            units = eye(4);
            for j = 1:4
                for k = 1:count
                    if strcmp(side,'left')
                        F.generators(:,k,j) = localQmul(units(:,j),P1(:,k));
                    else
                        F.generators(:,k,j) = localQmul(P1(:,k),units(:,j));
                    end
                end
            end
            F.bounds = localBounds(data.bounds,4);
            F.degree = count-1;
            F.parameterNames = {'real(t)';'i(t)';'j(t)';'k(t)'};
            F.side = side;

        otherwise
            error('qstab:familycomponents:UnknownUncertainty', ...
                'Unknown uncertainty kind "%s".',kind);
    end
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end

function C = localPolynomial(P)
    C = qstab.polycomponents(P);
    if ~isnumeric(C) || ~isreal(C) || size(C,1) ~= 4 || ...
            isempty(C) || any(~isfinite(C(:)))
        error('qstab:familycomponents:InvalidPolynomial', ...
            'A scalar polynomial must have finite real quaternion components.');
    end
    C = double(C);
end

function B = localBounds(value,m)
    if ~isnumeric(value) || ~isreal(value) || ...
            ~isequal(size(value),[m,2]) || any(~isfinite(value(:)))
        error('qstab:familycomponents:InvalidBounds', ...
            'bounds must be a finite real %d-by-2 array.',m);
    end
    B = double(value);
    if any(B(:,1)>B(:,2))
        error('qstab:familycomponents:BoundsOrder', ...
            'Each lower parameter bound must be <= its upper bound.');
    end
end

function localRequiredStruct(value,fields,kind)
    if ~isstruct(value) || ~isscalar(value)
        error('qstab:familycomponents:InvalidStruct', ...
            'The %s family must be a scalar struct.',kind);
    end
    for j = 1:numel(fields)
        if ~isfield(value,fields{j})
            error('qstab:familycomponents:MissingField', ...
                'The %s family needs field %s.',kind,fields{j});
        end
    end
end

function names = localParameterNames(m,prefix)
    names = cell(m,1);
    for j = 1:m
        names{j} = sprintf('%s(%d)',prefix,j);
    end
end

function tf = localIsText(value)
    tf = (ischar(value) && (isrow(value) || isempty(value))) || ...
         (isa(value,'string') && isscalar(value));
end

function c = localQmul(a,b)
    c = [a(1)*b(1)-a(2:4).'*b(2:4); ...
         a(1)*b(2:4)+b(1)*a(2:4)+cross(a(2:4),b(2:4))];
end
