function [interval,Info] = raymargin(X0, X1, u, varargin)
%RAYMARGIN  First loss of strict stability along one quaternion direction.
%
%   [ab,Info] = qstab.raymargin(F,u) uses a qstab.family('oneparameter',...)
%   definition. ab is always the FULL-law first-exit interval. The distinct
%   Info.parameterDomain.admissibleAlpha is the closed intersection of the
%   given box with t=alpha*u (empty if missed); its edges are not contacts.
%   ab = qstab.raymargin(P0,P1,u)
%   [ab,Info] = qstab.raymargin(P0,P1,u,'Side','right')
%   ab = qstab.raymargin(P0,P1,u,'Variable','z')
%   ab = qstab.raymargin(A0,A1,u,'Type','matrix','Domain','discrete')
%   Any form accepts the name/value options below.
%
%   Result. ab is a 1-by-2 real row [alphaMinus alphaPlus]. Starting at
%   alpha=0, ab(1) is the first contact in the negative direction and ab(2)
%   in the positive direction. If both are resolved, the nominal member's
%   connected interval of strict stability and regularity is (ab(1),ab(2)).
%   Exterior points may be stable again. The endpoints are signed alpha,
%   not magnitudes of the perturbation t=alpha*u. A contact at alpha=a has
%   |t|=abs(a)*abs(u). The direction u is not normalized by this function.
%   +/-Inf means no finite exit was found for that resolved direction;
%   NaN means the endpoint is unresolved. Read Info.status,
%   Info.negativeStatus, Info.positiveStatus and Info.message if needed.
%
%   Polynomial inputs. P0 is the strictly stable nominal scalar polynomial.
%   P1 is the perturbation polynomial. Supply coefficients in ASCENDING
%   powers: [c0 c1 ... cn] for numeric or native quaternion coefficients,
%   or a struct with real component arrays a,b,c,d of equal size. u is
%   one finite nonzero scalar quaternion, numeric scalar or component struct.
%   For continuous time the nominal degree is at least one and deg(P1)<deg(P0).
%   For discrete time the fixed declared grade is kept; loss of the forward
%   leading coefficient is a regularity event, distinct from a finite zero.
%
%   Matrix inputs. A0 and A1 are square quaternion or numeric matrices
%   of the same size. Select Type='matrix' explicitly. The nominal A0
%   must have strictly stable right spectrum. u has the same meaning.
%
%   Options (name followed by value):
%     Side          'left' (default): P0+t*P1 or A0+(t*I)*A1;
%                   'right': P0+P1*t or A0+A1*(t*I).
%     Type          'polynomial' (default) or 'matrix'.
%     Variable      's' (default), 'z' (forward shift), 'zinv' (delay).
%     Domain        'continuous' (default) or 'discrete' for state matrices.
%     Grade         [] infers the grade; discrete polynomials only.
%     Tolerance     1e-9; numerical accuracy/conditioning guard, 0<Tol<1.
%     VerifyBoundary true; check computed boundary contacts numerically.
%     MaxGuardCells 2047; discrete verification cell budget.
%     MaxGuardDepth 50; discrete verification subdivision depth limit.
%     MaxGuardTime  5 seconds; discrete verification time budget.
%
%   Example: first exit followed by re-entry, for a real direction u=1.
%     P0 = [0.5 1 1 1]; P1 = [4 1 1];
%     [ab,Info] = qstab.raymargin(P0,P1,1);
%     alphaMinus = ab(1); alphaPlus = ab(2);
%     % Expected: -1/8 and 1-1/sqrt(2). alpha=2 is stable again.
%     [tf,S] = qstab.isstable(P0 + 2*[P1 0]);  % tf=1
%     disp([alphaMinus alphaPlus]); disp(S.status);
%
%   Numerical details. Pole/zero classes are tested in ordinary floating
%   point. Ambiguous earlier contacts are not skipped in favour of a later
%   loss. The discrete dense ray method supports orders/grades up to 16.
%   See also qstab.freqmargin, qstab.stabilityradius, qstab.isstable.
try
if isstruct(X0)&&isscalar(X0)&&all(isfield(X0,{'p0','p1','bounds','side'}))&&~isfield(X0,'schema')
    % Legacy struct fields do not carry variable or grade; explicit options can set them.
    if nargin>=3,a=[{u} varargin];else,a=varargin;end
    meta={};
    if mod(numel(a),2),error('qstab:family:Options','Supply complete name/value pairs.');end
    for j=1:2:numel(a)
        if (ischar(a{j})||isstring(a{j}))&&any(strcmpi(a{j},{'Variable','Grade'})),meta(end+1:end+2)=a(j:j+1);end
    end
    X0=qstab.family('oneparameter',X0.p0,X0.p1,'Side',X0.side,'Bounds',X0.bounds,meta{:});
end
if qstab.internal.isFamilyDefinition(X0)
    if nargin<2,error('qstab:family:Direction','Supply a nonzero quaternion direction.');end
    if nargin>=3,varargin=[{u},varargin];end
    [D,opts,dom]=qstab.internal.familyMargin(X0,varargin,'raymargin');
    Info=qstab.internal.raymarginReport(D.p0,D.p1,X1,opts{:});
    [dom.admissibleAlpha,dom.nominalAdmissible,dom.alphaRangeStatus]= ...
        qstab.internal.parameterRay(D.bounds,Info.directionComponents);
    Info.parameterDomain=dom;
else
    Info = qstab.internal.raymarginReport(X0, X1, u, varargin{:});
end
interval = Info.interval;
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
