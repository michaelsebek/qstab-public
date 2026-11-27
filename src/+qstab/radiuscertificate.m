function [bounds,Info] = radiuscertificate(P0,P1,rBounds,varargin)
%RADIUSCERTIFICATE Verify rational lower/upper bounds for a quaternion radius.
%   [bounds,Info] = qstab.radiuscertificate(F,rBounds)
%   [bounds,Info] = qstab.radiuscertificate(P0,P1,rBounds)
%   [bounds,Info] = qstab.radiuscertificate(P0,P1,rBounds,'Side','right')
%
%   F is an unrestricted qstab.family('oneparameter',...) definition. P0 is
%   a Hurwitz scalar polynomial and degree(P1)<degree(P0). One constant
%   quaternion t forms P0+t*P1 for Side='left', or P0+P1*t for Side='right'.
%   Coefficient vectors are ascending. Native quaternion values, numeric
%   coefficients and ordinary component structs a,b,c,d are accepted.
%
%   rBounds is [lower upper], or a two-element cell array of exact rational
%   text, for example {'378937381963/1000000000000', ...
%   '378937381964/1000000000000'}. The lower endpoint may be zero. Otherwise
%   both endpoints must be finite, rational and 0<=lower<upper.
%
%   bounds is the verified 1-by-2 enclosure, or [NaN NaN] if a requested
%   proof fails. Inf Inf is returned for P1=0. Double endpoints are rounded
%   outwards; Info.exactBounds contains the authoritative rational strings.
%   Info.certified and Info.ok are true only when the nominal assumptions
%   and BOTH bounds were proved in this invocation. Failed certificates
%   are not instability decisions. No numerical search is silently used.
%
%   The proof covers ALL frequencies, using x=omega^2 and exact rational
%   Sturm root counts of Q_r=delta0-2*r^2*delta01+r^4*delta1. A positive lower
%   endpoint is accepted if Q_lower has no nonnegative real root. For the
%   upper endpoint, a root of Q_upper proves that a critical gain of norm
%   no greater than upper exists. Alternatively UpperX supplies one exact
%   nonnegative rational x at which Q_upper<=0 or T_upper=n0-upper^2*n1<=0.
%   Repeated roots and flat extrema are included by square-free reduction.
%
%   Options:
%     Side            'left' (default) or 'right'.
%     RationalData    [] (default), or struct('P0',C0,'P1',C1), where C0,C1
%                     are 4-by-N ascending real rational component matrices.
%                     Entries may be numeric, symbolic rational, or fraction
%                     strings in cells. Explicit data must round to the
%                     supplied input coefficients, without shape changes.
%     UpperX          []: count roots of Q_upper. Otherwise certify the
%                     supplied x=omega^2 directly; it is NOT omega.
%     Variable        's' only. Grade must be [].
%     ExtendParameter false. A finite family domain requires explicit true;
%                     the result then concerns its law on all quaternion t.
%
%   Exact numeric semantics: a supplied double denotes its exact binary
%   value (sym(value,'f')); no rational approximation is guessed. Use
%   fraction text in RationalData when the intended source is, for example,
%   exactly 1/10 rather than its stored binary approximation. Symbolic input
%   coefficients must be rational. Returned reports contain no sym objects.
%
%   Example (requires Symbolic Math Toolbox):
%     [b,I] = qstab.radiuscertificate([1 1 1],1,{'0.8660','0.8661'});
%     assert(I.certified); b
%
%   Requirement: optional Symbolic Math Toolbox. Ordinary numerical calls
%   to qstab.stabilityradius/freqmargin remain independent of that toolbox.
%   This function does not modify files, paths, symbolic assumptions or RNG.
%   See also qstab.raycertificate, qstab.stabilityradius, qstab.freqmargin.
try
    if qstab.internal.isFamilyDefinition(P0)
        if nargin<2,error('qstab:exact:Arguments','Supply the requested radius bounds.');end
        if nargin>=3,varargin=[{rBounds},varargin];end
        [D,args,domain]=qstab.internal.familyMargin(P0,varargin,'radiuscertificate');
        [bounds,Info]=qstab.radiuscertificate(D.p0,D.p1,P1,args{:});
        Info.parameterDomain=domain;return
    end
    if nargin<3,error('qstab:exact:Arguments','Supply P0, P1 and two candidate radius bounds.');end
    o=options(varargin{:});qstab.internal.exactOneParameter('require');
    [A,B,~,source]=qstab.internal.exactOneParameter('data',P0,P1,[],o.RationalData);
    r=qstab.internal.exactOneParameter('rational',rBounds);r=reshape(r,1,numel(r));
    if numel(r)~=2||signExact(r(1))<0||signExact(r(2)-r(1))<=0
        error('qstab:exact:Bounds','Use two finite rational endpoints with 0<=lower<upper.');
    end
    bounds=[NaN NaN];Info=baseInfo(o.Side,source);
    Info.requestedBounds=qstab.internal.exactOneParameter('strings',r);
    Info.inputComponents=struct('P0',{qstab.internal.exactOneParameter('strings',A)}, ...
        'P1',{qstab.internal.exactOneParameter('strings',B)});
    if size(A,2)<2||(~allZero(B)&&size(B,2)>=size(A,2))
        Info.status='assumptionFailed';Info.message='Require degree(P0)>=1 and degree(P1)<degree(P0).';return
    end
    [nominalOK,Info.nominal]=qstab.internal.exactOneParameter('nominal',A);
    if ~nominalOK
        Info.status='assumptionFailed';Info.message='The nominal norm polynomial is not strictly Hurwitz over the exact input field.';return
    end
    if allZero(B)
        bounds=[Inf Inf];Info.bounds=bounds;Info.exactBounds={'Inf','Inf'};
        Info.ok=true;Info.certified=true;Info.validated=true;Info.status='unbounded';
        Info.message='The perturbation polynomial is exactly zero.';return
    end
    [n0,n1,d0,d1,d01,x]=qstab.internal.exactOneParameter('frequency',A,B,o.Side);
    Qlo=expand(d0-2*r(1)^2*d01+r(1)^4*d1);
    Qhi=expand(d0-2*r(2)^2*d01+r(2)^4*d1);Thi=expand(n0-r(2)^2*n1);
    Info.polynomials=struct('variable','x=omega^2','symbol',char(x), ...
        'delta0',char(d0),'delta1',char(d1),'delta01',char(d01), ...
        'Qlower',char(Qlo),'Qupper',char(Qhi),'Tupper',char(Thi));
    if signExact(r(1))==0
        Info.lowerProof=struct('accepted',true,'method','nonnegativity-of-radius', ...
            'nonnegativeRootCount',NaN);
    else
        L=countNonnegative(Qlo,x);
        L.accepted=L.nonnegativeRootCount==0;Info.lowerProof=L;
    end
    if isempty(o.UpperX)
        U=countNonnegative(Qhi,x);U.accepted=U.nonnegativeRootCount>0;
        Info.upperProof=U;
    else
        ux=qstab.internal.exactOneParameter('rational',o.UpperX);
        if ~isscalar(ux)||signExact(ux)<0,error('qstab:exact:UpperX','UpperX must be one nonnegative rational x=omega^2.');end
        qx=expand(subs(Qhi,x,ux));tx=expand(subs(Thi,x,ux));
        Info.upperProof=struct('accepted',signExact(qx)<=0||signExact(tx)<=0, ...
            'method','exact-rational-frequency-witness','x',char(ux), ...
            'Qvalue',char(qx),'Tvalue',char(tx));
    end
    Info.certified=Info.lowerProof.accepted&&Info.upperProof.accepted;
    Info.ok=Info.certified;Info.validated=Info.certified;
    if Info.certified
        bounds=qstab.internal.exactOneParameter('enclosure',r(1),r(2));
        Info.bounds=bounds;Info.exactBounds=qstab.internal.exactOneParameter('strings',r);
        Info.status='certified';Info.message='Both rational endpoints are proved over the complete frequency half-line.';
    else
        Info.status='notCertified';Info.message='At least one supplied bound failed its exact certificate; inspect lowerProof and upperProof.';
    end
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'),throwAsCaller(problem);end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end

function P=countNonnegative(q,x)
atZero=signExact(subs(q,x,0))==0;
[S,f]=qstab.internal.exactOneParameter('sturm',q,x);
[n,v0,vi]=qstab.internal.exactOneParameter('count',S,x,sym(0),Inf);
P=struct('accepted',false,'method','exact-square-free-Sturm-root-count', ...
    'positiveRootCount',n,'rootAtZero',atZero,'nonnegativeRootCount',n+double(atZero), ...
    'variationAtZero',v0,'variationAtInfinity',vi,'squareFreePolynomial',char(f), ...
    'sturmSequence',{cellfun(@char,S,'UniformOutput',false)});
end

function I=baseInfo(side,source)
I=struct('ok',false,'certified',false,'validated',false,'status','undetermined', ...
    'assurance','exact-rational','computedThisRun',true,'method','rational-radius-polynomial-certificate', ...
    'dependency','Symbolic Math Toolbox','inputSemantics',source,'side',side,'variable','s', ...
    'bounds',[NaN NaN],'exactBounds',{{}},'nominal',struct(), ...
    'lowerProof',struct(),'upperProof',struct(),'message','', ...
    'parameterDomain',struct('scope','unrestrictedModel','isBounded',false));
end
function o=options(varargin)
o=struct('Side','left','RationalData',[],'UpperX',[],'Variable','s','Grade',[],'ExtendParameter',false);
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
if ~(isscalar(o.ExtendParameter)&&(islogical(o.ExtendParameter)||(isnumeric(o.ExtendParameter)&&isreal(o.ExtendParameter)&&any(o.ExtendParameter==[0 1]))))
    error('qstab:exact:Options','ExtendParameter must be a logical scalar.');
end
end
function s=signExact(v),s=qstab.internal.exactOneParameter('sign',v);end
function tf=allZero(a),tf=qstab.internal.exactOneParameter('iszero',a);end
