function D = randomFamilyCoefficients(o)
%RANDOMFAMILYCOEFFICIENTS Internal generator of 4-by-(g+1) real coefficient blocks.
% Columns are ascending in the requested variable. Rows are 1,i,j,k.
% Draw first in the forward variable; reverse EVERY block at the declared
% grade for delay input, including the perturbation's zero leading column.
% The caller validates the options and provides the selected random stream.
g = o.degree;
C = nominal(g,o);
switch o.kind
    case 'interval'
        % Centre plus/minus independently sampled, strictly positive
        % half-widths on every free coefficient component.
        mask = freeComponents(g);
        if o.integer
            r = o.width*o.drawint(2,nnz(mask),1);
        else
            r = o.width*(0.5+0.5*o.uniform(nnz(mask),1));
        end
        W = zeros(size(C)); W(mask) = r;
        L = C-W; U = C+W;
        if ~all(isfinite(L(:))) || ~all(isfinite(U(:))) || ...
                any(L(mask)>=C(mask)) || any(U(mask)<=C(mask))
            rangeError();
        end
        % Interval import checks endpoint lengths as well as Grade. If an
        % entire last delay block cancelled exactly, use a symmetric
        % expansion of one component interval to keep both endpoint lengths.
        % With g>=1 a single component suffices. Enlarge both sides equally.
        if strcmp(o.variable,'zinv') && g>0 && (~any(L(:,1)) || ~any(U(:,1)))
            ix = find(W(:,1)>0,1);
            W(ix,1) = 2*W(ix,1);
            L = C-W; U = C+W;
            if ~all(isfinite(L(:))) || ~all(isfinite(U(:))), rangeError(); end
        end
        D = struct('lower',orient(L,o.variable),'upper',orient(U,o.variable));
    case 'polytope'
        V = cell(1,o.count);
        V{1} = orient(C,o.variable);
        for n=2:o.count
            A = increment(g,o);
            Vn = C+A;
            if ~all(isfinite(Vn(:))) || ~any(Vn(:)~=C(:)), rangeError(); end
            V{n} = orient(Vn,o.variable);
        end
        D = struct('vertices',{V});
    case 'affine'
        G = cell(1,o.count);
        for n=1:o.count
            G{n} = orient(increment(g,o),o.variable);
        end
        D = struct('center',orient(C,o.variable),'generators',{G}, ...
            'bounds',o.bounds);
    case 'oneparameter'
        D = struct('p0',orient(C,o.variable), ...
            'p1',orient(increment(g,o),o.variable), ...
            'bounds',o.bounds,'side',o.side);
    otherwise
        error('qstab:randfamily:Kind','Unknown family kind.');
end
end

function C = nominal(g,o)
v=o.variable; intMode=o.integer;
if g==0
    % Every generated constant family lies in Re(P)=1.
    if intMode, b = o.drawint(5,3,1)-3; else, b = 2*o.uniform(3,1)-1; end
    C = [1;b];
elseif strcmp(v,'s')
    % Central polynomial indeterminate. Multiply on the LEFT by s+a at
    % each step using Hamilton multiplication. Re(a)>0; the algebraic
    % norm is the product of the real quadratic norms of these factors.
    C = [1;0;0;0];
    for n=1:g
        if intMode
            a = [o.drawint(2);o.drawint(3,3,1)-2];
        else
            a = [0.75+0.5*o.uniform(); 2*o.uniform(3,1)-1];
        end
        ar=a(1); bi=a(2); cj=a(3); dk=a(4);
        L = [ar -bi -cj -dk; bi ar -dk cj; cj dk ar -bi; dk -cj bi ar];
        C = [L*C,zeros(4,1)] + [zeros(4,1),C];
        if ~all(isfinite(C(:))), rangeError(); end
    end
else
    C = zeros(4,g+1); C(1,end) = 1;
    if ~intMode
        % For |lambda|>=1, the sum of lower-order terms is strictly
        % smaller in norm than the monic term (in the complex adjoint too).
        A = 2*o.uniform(4,g)-1;
        if ~any(A(:)), A(2,1)=1; end
        C(:,1:g) = ((0.2+0.4*o.uniform())/sum(sqrt(sum(A.^2,1))))*A;
    end
    % In integer mode the Schur nominal is z^g (zero forward modes).
end
end

function mask = freeComponents(g)
mask = true(4,g+1);
if g==0, mask(1)=false; else, mask(:,end)=false; end
end

function A = increment(g,o)
w=o.width; intMode=o.integer;
mask = freeComponents(g);
A = zeros(size(mask));
if intMode
    A(mask) = w*(o.drawint(5,nnz(mask),1)-3);
else
    A(mask) = w*(2*o.uniform(nnz(mask),1)-1);
end
% Nonzero direction even in integer mode. The highest *free* block has
% a nonreal component, allowing tests that distinguish left and right.
j = max(1,g);
if ~any(A(2:4,j)), A(2,j)=w; end
if ~all(isfinite(A(:))), rangeError(); end
end

function B = orient(A,v)
if strcmp(v,'zinv'), B=A(:,end:-1:1); else, B=A; end
end

function rangeError()
error('qstab:randfamily:CoefficientRange', ...
    'Degree or Width is outside the representable coefficient range for this draw.');
end
