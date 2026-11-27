function [gammaBounds,Info] = raycertificate(P0,P1,d,varargin)
%RAYCERTIFICATE Prove the first positive loss on one quaternion ray exactly.
%   [gammaBounds,Info] = qstab.raycertificate(F,d)
%   [gammaBounds,Info] = qstab.raycertificate(P0,P1,d)
%   [gammaBounds,Info] = qstab.raycertificate(P0,P1,d,'Side','right')
%
%   F is an unrestricted qstab.family('oneparameter',...) definition. The
%   raw form describes P0+t*P1 (Side='left', default) or P0+P1*t ('right'),
%   with t=gamma*d and gamma>0. d is one NONZERO quaternion and is never
%   normalized. P0 must be Hurwitz, degree(P0)>=1 and degree(P1)<degree(P0).
%   Coefficients are ascending numeric/native quaternion vectors or ordinary
%   structs of the four real components a,b,c,d. Only continuous time is used.
%
%   gammaBounds=[lo hi] encloses the FIRST POSITIVE real gamma where strict
%   stability is lost, including isolated tangential contacts. It is not
%   the signed [negative positive] interval returned by qstab.raymargin.
%   Call with -d for the opposite direction. [Inf Inf] proves no positive
%   contact. A finite root outside the isolation budget gives [NaN NaN].
%   Exact rational endpoints are Info.exactBounds; returned doubles are
%   rounded outwards. The gain norm is gamma*abs(d), not gamma itself.
%
%   The real norm polynomial is quadratic in gamma. Its reduced Hurwitz
%   determinant is formed exactly, then its square-free part is isolated
%   with exact Sturm counts. Even-multiplicity roots are retained. No
%   floating-point eigenvalue is accepted as a proof or rounded to real.
%   The nominal Hurwitz property and fixed degree are also checked exactly.
%
%   Options:
%     Side            'left' (default) or 'right'.
%     RationalData    [] or struct('P0',C0,'P1',C1,'Direction',d4).
%                     C0,C1 are 4-by-N ascending rational component arrays;
%                     d4 is 4-by-1. Numeric, sym rational and cell fraction
%                     strings are accepted. Direction may be omitted.
%                     Explicit data must round to the supplied coefficients.
%     AbsTol          1e-12, requested absolute gamma enclosure width.
%     MaxBisections   256, positive integer. Exact counting is independent
%                     of this precision budget; inspect widthTargetMet.
%     VerifyLifting   true: verify the contact on the original polynomial
%                     where an exact algebraic representation is available.
%                     Lifting status is separate from the exact first-root
%                     proof and is never replaced by a numerical residual.
%     Variable        's' only. Grade must be [].
%     ExtendParameter false; finite family bounds require explicit true.
%
%   Info.certified/ok is true only after an exact first-root proof (or an
%   exact proof of no positive root). Info.isolation stores the determinant,
%   square-free polynomial, Sturm sequence and root counts. Reports contain
%   ordinary MATLAB data and rational strings, so they can be saved normally.
%   Numeric doubles are interpreted as exact stored binary rationals using
%   sym(value,'f'); decimal/fraction source data belong in RationalData.
%
%   Example (a tangential contact, requires Symbolic Math Toolbox):
%     P0=struct('a',[1 2 1],'b',[0 0 0],'c',[1 0 0],'d',[0 0 0]);
%     P1=struct('a',[-1 -1],'b',[0 0],'c',[-1 -1],'d',[0 0]);
%     [ab,I]=qstab.raycertificate(P0,P1,1); assert(I.certified); ab
%     % The exact first positive contact is gamma=1, despite re-entry.
%
%   Optional dependency: Symbolic Math Toolbox. Numerical raymargin remains
%   available without it. Symbolic determinant computations may be costly.
%   See also qstab.raymargin, qstab.radiuscertificate, qstab.freqmargin.
try
    if qstab.internal.isFamilyDefinition(P0)
        if nargin<2,error('qstab:exact:Arguments','Supply a quaternion ray direction.');end
        if nargin>=3,varargin=[{d},varargin];end
        [D,args,domain]=qstab.internal.familyMargin(P0,varargin,'raycertificate');
        [gammaBounds,Info]=qstab.raycertificate(D.p0,D.p1,P1,args{:});
        Info.parameterDomain=domain;return
    end
    if nargin<3,error('qstab:exact:Arguments','Supply P0, P1 and one quaternion direction.');end
    o=options(varargin{:});qstab.internal.exactOneParameter('require');
    [A,B,d4,source]=qstab.internal.exactOneParameter('data',P0,P1,d,o.RationalData);
    gammaBounds=[NaN NaN];Info=baseInfo(o.Side,source);
    Info.inputComponents=struct('P0',{qstab.internal.exactOneParameter('strings',A)}, ...
        'P1',{qstab.internal.exactOneParameter('strings',B)}, ...
        'Direction',{qstab.internal.exactOneParameter('strings',d4)});
    Info.directionComponents=double(d4);Info.directionNormSquared=char(sum(d4.^2));
    if size(A,2)<2||(~allZero(B)&&size(B,2)>=size(A,2))
        Info.status='assumptionFailed';Info.message='Require degree(P0)>=1 and degree(P1)<degree(P0).';return
    end
    [nominalOK,Info.nominal]=qstab.internal.exactOneParameter('nominal',A);
    if ~nominalOK
        Info.status='assumptionFailed';Info.message='The nominal norm polynomial is not strictly Hurwitz over the exact input field.';return
    end
    [h,N,g,s,C]=qstab.internal.exactOneParameter('ray',A,B,d4,o.Side);
    Info.determinant=char(h);Info.determinantVariable=char(g);Info.determinantCoefficients= ...
        qstab.internal.exactOneParameter('strings',fliplr(coeffs(h,g,'All')));
    [lo,hi,Info.isolation]=qstab.internal.exactOneParameter('isolate',h,g,o.AbsTol,o.MaxBisections);
    Info.candidateBounds=qstab.internal.exactOneParameter('enclosure',lo,hi);
    if ~Info.isolation.isolated
        Info.message='The budget did not isolate the first positive root from its neighbours.';return
    end
    gammaBounds=Info.candidateBounds;Info.bounds=gammaBounds;
    Info.ok=true;Info.certified=true;Info.validated=true;
    Info.widthTargetMet=Info.isolation.widthTargetMet;
    if isnumeric(lo)&&isinf(lo)
        Info.status='unbounded';Info.exactBounds={'Inf','Inf'};
        Info.message='Exact Sturm counting proves that no positive first contact exists.';return
    end
    Info.status='certified';Info.exactBounds={char(lo),char(hi)};
    Info.widthTargetMet=Info.isolation.widthTargetMet;
    Info.message='The first positive root is isolated exactly, including all multiple contacts.';
    if o.VerifyLifting
        Info.lifting=qstab.internal.exactOneParameterLifting(C,N,g,s,h,lo,hi);
        if strcmp(Info.lifting.status,'inconsistent')
            error('qstab:exact:Internal','The exact first-root proof and original-polynomial lifting are inconsistent.');
        end
    end
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'),throwAsCaller(problem);end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end

function I=baseInfo(side,source)
I=struct('ok',false,'certified',false,'validated',false,'status','undetermined', ...
    'assurance','exact-rational','computedThisRun',true,'method','exact-norm-Hurwitz-first-root', ...
    'dependency','Symbolic Math Toolbox','inputSemantics',source,'side',side,'variable','s', ...
    'bounds',[NaN NaN],'exactBounds',{{}},'nominal',struct(),'isolation',struct(), ...
    'widthTargetMet',false,'message','', ...
    'lifting',struct('status','notRequested','exact',false,'numericalWitness',false), ...
    'parameterDomain',struct('scope','unrestrictedModel','isBounded',false));
end
function o=options(varargin)
o=struct('Side','left','RationalData',[],'AbsTol',1e-12,'MaxBisections',256, ...
    'VerifyLifting',true,'Variable','s','Grade',[],'ExtendParameter',false);
if mod(numel(varargin),2),error('qstab:exact:Options','Use complete name/value pairs.');end
names=fieldnames(o);
for k=1:2:numel(varargin)
    name=varargin{k};if ~(ischar(name)||(isstring(name)&&isscalar(name))),error('qstab:exact:Options','Option names must be text.');end
    j=find(strcmpi(char(name),names),1);if isempty(j),error('qstab:exact:Options','Unknown option %s.',char(name));end
    o.(names{j})=varargin{k+1};
end
if ~((ischar(o.Side)&&isrow(o.Side))||(isstring(o.Side)&&isscalar(o.Side)))
    error('qstab:exact:Side','Side must be left or right.');
end
if ~((ischar(o.Variable)&&isrow(o.Variable))||(isstring(o.Variable)&&isscalar(o.Variable)))
    error('qstab:exact:Variable','Variable must be s.');
end
o.Side=lower(char(o.Side));
if ~ismember(o.Side,{'left','right'}),error('qstab:exact:Side','Side must be left or right.');end
if ~strcmpi(o.Variable,'s')||~isempty(o.Grade),error('qstab:exact:Variable','Exact one-parameter certificates currently require Variable=s and Grade=[].');end
if ~(isnumeric(o.AbsTol)&&isreal(o.AbsTol)&&isscalar(o.AbsTol)&&isfinite(o.AbsTol)&&o.AbsTol>0)
    error('qstab:exact:Options','AbsTol must be a positive finite scalar.');
end
if ~(isnumeric(o.MaxBisections)&&isreal(o.MaxBisections)&&isscalar(o.MaxBisections)&& ...
        isfinite(o.MaxBisections)&&o.MaxBisections>=1&&o.MaxBisections==floor(o.MaxBisections))
    error('qstab:exact:Options','MaxBisections must be a positive integer.');
end
for name={'VerifyLifting','ExtendParameter'}
    v=o.(name{1});
    if ~(isscalar(v)&&(islogical(v)||(isnumeric(v)&&isreal(v)&&any(v==[0 1]))))
        error('qstab:exact:Options','%s must be a logical scalar.',name{1});
    end
end
end
function tf=allZero(a),tf=qstab.internal.exactOneParameter('iszero',a);end
