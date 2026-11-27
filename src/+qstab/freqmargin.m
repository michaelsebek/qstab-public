function [rho, info] = freqmargin(P0, P1, omega, varargin)
%FREQMARGIN  Smallest |t| producing a boundary zero at a specified frequency.
%   [rho,I]=qstab.freqmargin(F,omega,...) accepts a oneparameter qstab.family.
%   These are all-H margins. For finite Bounds use ExtendParameter=true
%   only if you intend to extend the law beyond the specified domain.
%
%   rho = qstab.freqmargin(P0,P1,omega)
%   [rho,Info] = qstab.freqmargin(P0,P1,omega,'Side','right')
%   [rho,Info] = qstab.freqmargin(P0,P1,theta,'Variable','z')
%
%   P0 is a nominal strictly stable scalar polynomial. P1 is its perturbation
%   polynomial. Enter numeric/native quaternion coefficients in ascending
%   order or use a struct of real component arrays a,b,c,d. In continuous
%   time require deg(P0)>=1 and deg(P1)<deg(P0). One constant quaternion t
%   multiplies P1 on the chosen side for all coefficients.
%   omega is a nonnegative real scalar or vector; each query tests the
%   whole imaginary quaternion sphere of radius omega(k). In discrete time
%   supply angles theta in [0,pi] radians in the forward z variable.
%
%   rho has the shape of the query. rho(k) is the minimum parameter norm
%   for that class, not a global minimum over all frequencies. Inf denotes
%   an unattainable contact when the calculation establishes this; check
%   Info.valid(k) before interpreting a result. For a finite witness,
%   Info.parameterComponents(:,k) and Info.boundaryComponents(:,k) are
%   four real components in the order real,i,j,k.
%
%   Options:
%     Side         'left' (default): P0+t*P1; 'right': P0+P1*t.
%     Variable     's' (default), 'z' (forward), 'zinv' (delay).
%     Grade        [] infers the declared grade; discrete polynomials only.
%     Tolerance    1e-10, positive and less than one; numerical checks.
%     FrequencyScale 'auto' (default), or a positive real tau. For s only,
%                  compute with s=tau*sigma and a common positive coefficient
%                  scale. omega and every returned boundary point remain in
%                  the ORIGINAL s units; t and rho are unchanged. Auto uses
%                  a power-of-two tau when coefficient units are unbalanced.
%     CheckNominal true; false skips the nominal stability check when its
%                  result is already known.
%   Info.scaling records this change of units. Info.geometry{k} describes
%   the nearest gains at one frequency: 'point', 'sphere', 'undetermined',
%   or 'empty'. A sphere has centerComponents (4-by-1), basisComponents
%   (4-by-3, orthonormal) and sphereRadius, so its gains are center+r*basis*v
%   for real norm(v)=1. A small numerical gap alone NEVER establishes a
%   sphere: that label requires the structural common-orientation test.
%   Info.svdRadius and Info.svdParameterComponents provide an independent
%   ordered 2-by-2 complex-adjoint comparison, using base MATLAB only.
%   Info.zeroWitnesses{k} gives both evaluation-side zeros when the central
%   quadratic remainder permits reliable recovery. These numerical details
%   are not exact multiplicity or whole-frequency certificates.
%
%   Example:
%     [rho,I] = qstab.freqmargin([2 1],1,[0 1 2]);
%     disp(rho);  % sqrt(4 + [0 1 2].^2)
%     t0 = I.parameterComponents(:,1);
%
%   See also qstab.stabilityradius, qstab.raymargin.

try
    if qstab.internal.isFamilyDefinition(P0)
        if nargin<2,error('qstab:family:Frequency','Supply the requested frequency.');end
        if nargin>=3,varargin=[{omega},varargin];end
        [D,opts,dom]=qstab.internal.familyMargin(P0,varargin,'freqmargin');
        [rho,info]=qstab.freqmargin(D.p0,D.p1,P1,opts{:});info.parameterDomain=dom;return
    end
    o=struct('Side','left','Tolerance',1e-10,'CheckNominal',true,'Variable','s','Grade',[], ...
        'FrequencyScale','auto');
    o=parseopts(o,varargin{:}); side=lower(char(o.Side));
    if ~any(strcmp(side,{'left','right'}))
        error('qstab:freqmargin:Side','Side must be left or right.');
    end
    if ~(isnumeric(omega)&&isreal(omega)&&isvector(omega)&&~isempty(omega)&&...
            all(isfinite(omega(:)))&&all(omega(:)>=0))
        error('qstab:freqmargin:Frequency','omega must be a nonempty finite nonnegative real vector.');
    end
    if ~(isnumeric(o.Tolerance)&&isreal(o.Tolerance)&&isscalar(o.Tolerance)&&...
            isfinite(o.Tolerance)&&o.Tolerance>0&&o.Tolerance<1)
        error('qstab:freqmargin:Tolerance','Tolerance must lie strictly between 0 and 1.');
    end
    if ~(isscalar(o.CheckNominal)&&(islogical(o.CheckNominal)||...
            (isnumeric(o.CheckNominal)&&any(o.CheckNominal==[0,1]))))
        error('qstab:freqmargin:CheckNominal','CheckNominal must be a logical scalar.');
    end
    if any(strcmpi(o.Variable,{'z','zinv'}))
        if ~isAuto(o.FrequencyScale)
            error('qstab:freqmargin:FrequencyScale','FrequencyScale is available only for Variable=''s''.');
        end
        [rho,info]=qstab.internal.unitfrequency(P0,P1,omega,o); return
    elseif ~strcmpi(o.Variable,'s')
        error('qstab:freqmargin:Variable','Variable must be s, z or zinv.');
    end
    if ~isempty(o.Grade)
        error('qstab:freqmargin:Grade','Grade is currently a discrete-frequency option only.');
    end
    C0=qstab.polycomponents(P0); C1=qstab.polycomponents(P1);
    n0=find(any(C0~=0,1),1,'last'); n1=find(any(C1~=0,1),1,'last');
    if isempty(n0)||n0<2||(~isempty(n1)&&n1>=n0)
        error('qstab:freqmargin:Degree','Require degree(P0)>=1 and degree(P1)<degree(P0).');
    end
    C0=C0(:,1:n0);
    if isempty(n1), C1=zeros(4,1); else, C1=C1(:,1:n1); end
    [C0,C1,scaling]=qstab.internal.oneParameterScaling(C0,C1,o.FrequencyScale);
    workingOmega=double(omega)/scaling.frequencyScale;
    if any(~isfinite(workingOmega(:)))||any((omega(:)~=0)&(workingOmega(:)==0))
        error('qstab:freqmargin:FrequencyScale','A requested frequency is outside the representable range after scaling.');
    end
    if o.CheckNominal
        nominal=qstab.internal.isstableReport(aspoly(C0),'Variable','s','Tolerance',o.Tolerance);
        if ~strcmp(nominal.status,'stable')
            error('qstab:freqmargin:Nominal','P0 must be numerically verified Hurwitz.');
        end
    else
        nominal=struct('status','notChecked');
    end
    nw=numel(omega); rho=NaN(size(omega));
    info=struct('side',side,'variable','s','assurance','floating-point',...
        'valid',false(size(omega)),'parameterComponents',NaN(4,nw),...
        'boundaryComponents',NaN(4,nw),'evaluationResidual',NaN(size(omega)),...
        'objectiveRadius',NaN(size(omega)),'centerLowerEstimate',NaN(size(omega)),...
        'formulaWitnessRelativeDifference',NaN(size(omega)),...
        'minNominalModulus',NaN(size(omega)),'maxPerturbationModulus',NaN(size(omega)),...
        'nominal',nominal,'message',{cell(size(omega))}, ...
        'frequency',double(omega),'workingFrequency',workingOmega,'frequencyScale',scaling.frequencyScale, ...
        'scaling',scaling,'geometry',{cell(size(omega))},'zeroWitnesses',{cell(size(omega))}, ...
        'svdRadius',NaN(size(omega)),'svdRelativeDifference',NaN(size(omega)), ...
        'svdParameterComponents',NaN(4,nw),'singularValueGapRatio',NaN(size(omega)), ...
        'singularGapUncertaintyEstimate',NaN(size(omega)));
    [q0,common0]=commonorientation(C0);[q1,common1]=commonorientation(C1);
    for k=1:nw
        w=workingOmega(k);
        info.geometry{k}=emptygeometry('undetermined','The frequency calculation has not produced a verified witness.');
        info.zeroWitnesses{k}=emptyzeros();
        [x0,y0]=evenodd(C0,w); [x1,y1]=evenodd(C1,w);
        s0=max(abs([x0;y0])); s1=max(abs([x1;y1]));
        if ~isfinite(s0)||~isfinite(s1)||s0==0
            info.message{k}='Evaluation overflow or zero nominal evaluation.'; continue
        end
        x0=x0/s0; y0=y0/s0;
        a=x0'*x0+y0'*y0;
        if strcmp(side,'left'), c0=mulq(conjq(y0),x0); else, c0=mulq(x0,conjq(y0)); end
        c0=c0(2:4); m=norm(c0); gap=a-2*m;
        if gap<=o.Tolerance*max(a,realmin)
            info.message{k}='Nominal boundary separation is too small for this tolerance.'; continue
        end
        info.minNominalModulus(k)=s0*sqrt(gap);
        if ~isfinite(info.minNominalModulus(k))||info.minNominalModulus(k)==0
            info.message{k}='Nominal separation is outside the representable range.'; continue
        end
        if s1==0
            rho(k)=Inf; info.valid(k)=true; info.maxPerturbationModulus(k)=0;
            info.objectiveRadius(k)=Inf; info.centerLowerEstimate(k)=Inf;
            info.geometry{k}=emptygeometry('empty','The perturbation vanishes on the entire frequency sphere.');
            info.message{k}='Perturbation vanishes on this entire frequency sphere.'; continue
        end
        x1=x1/s1; y1=y1/s1;
        b=x1'*x1+y1'*y1;
        if strcmp(side,'left'), c1=mulq(conjq(y1),x1); else, c1=mulq(x1,conjq(y1)); end
        c1=c1(2:4);
        delta0=gap*(a+2*m); delta01=a*b-4*(c0'*c1);
        if m>0
            beta=(c0'*c1)/m; perp=c1-(beta/m)*c0;
            rad=4*((a*beta-b*m)^2+delta0*(perp'*perp));
        else
            rad=4*a*a*(c1'*c1);
        end
        den=delta01+sqrt(rad);
        if ~(isfinite(den)&&den>0)
            info.message{k}='Ill-conditioned quotient in the frequency formula.'; continue
        end
        scaleRatio=s0/s1;
        if ~isfinite(scaleRatio)||scaleRatio<realmin
            info.message{k}='Coefficient-scale ratio overflows or underflows working precision.'; continue
        end
        r2=delta0/den; r=sqrt(r2)*scaleRatio;
        v=-(c0-r2*c1); nv=norm(v);
        if nv<=100*eps*(norm(c0)+r2*norm(c1)+realmin), v=[1;0;0]; else, v=v/nv; end
        vq=[0;v];
        if strcmp(side,'left')
            e0=x0+mulq(y0,vq); e1=x1+mulq(y1,vq);
        else
            e0=x0+mulq(vq,y0); e1=x1+mulq(vq,y1);
        end
        en=e1'*e1;
        if en<=eps^2*max(b,realmin)||~isfinite(r)||r<realmin
            info.message{k}='Critical-parameter reconstruction is numerically singular.'; continue
        end
        if strcmp(side,'left')
            ts=-mulq(e0,conjq(e1))/en;
            err=e0+mulq(ts,e1);
        else
            ts=-mulq(conjq(e1),e0)/en;
            err=e0+mulq(e1,ts);
        end
        t=scaleRatio*ts;
        tsReturned=t/scaleRatio;
        if strcmp(side,'left'), err=e0+mulq(tsReturned,e1); else, err=e0+mulq(e1,tsReturned); end
        residual=norm(err)/max(realmin,norm(e0)+norm(ts)*norm(e1));
        disagreement=abs(norm(t)-r)/max(realmin,r);
        info.parameterComponents(:,k)=t;
        info.boundaryComponents(:,k)=double(omega(k))*vq;
        info.evaluationResidual(k)=residual;
        info.maxPerturbationModulus(k)=s1*sqrt(b+2*norm(c1));
        if residual>100*o.Tolerance||disagreement>100*o.Tolerance
            info.message{k}='Frequency formula and reconstructed critical parameter disagree.'; continue
        end
        rho(k)=norm(t); info.valid(k)=true;
        info.objectiveRadius(k)=r;
        info.formulaWitnessRelativeDifference(k)=disagreement;
        info.centerLowerEstimate(k)=max(0,min(r,rho(k))-abs(r-rho(k))-100*o.Tolerance*r);
        info.message{k}='Numerical nearest critical parameter at this frequency, with reconstructed witness.';
        [svdr,svdt,gapSvd,gapGuard]=svdcomparison(x0,y0,x1,y1,scaleRatio,side);
        info.svdRadius(k)=svdr;info.svdParameterComponents(:,k)=svdt;
        info.svdRelativeDifference(k)=abs(svdr-rho(k))/max(realmin,rho(k));
        info.singularValueGapRatio(k)=gapSvd;
        info.singularGapUncertaintyEstimate(k)=gapGuard;
        info.geometry{k}=gaingeometry(t,x0,y0,x1,y1,q0,q1,common0&&common1, ...
            scaleRatio,side,double(omega(k)),gapSvd,gapGuard,o.Tolerance);
        info.zeroWitnesses{k}=zerowitnesses(C0,C1,t,w,double(omega(k)),side,o.Tolerance);
    end
    % Preserve the public interpretation of these moduli in the caller's
    % coefficient units. They are diagnostics, not the scaled working data.
    info.minNominalModulus=pow2(info.minNominalModulus,scaling.coefficientScalePowerOfTwo);
    info.maxPerturbationModulus=pow2(info.maxPerturbationModulus,scaling.coefficientScalePowerOfTwo);
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end

function [x,y]=evenodd(C,w)
    x=zeros(4,1); y=zeros(4,1); z=-w*w;
    for k=size(C,2):-1:1
        if mod(k-1,2)==0, x=x*z+C(:,k); else, y=y*z+C(:,k); end
    end
    y=w*y;
end
function z=mulq(x,y)
    z=[x(1)*y(1)-x(2:4)'*y(2:4);...
        x(1)*y(2:4)+y(1)*x(2:4)+cross(x(2:4),y(2:4))];
end
function z=conjq(x)
    z=[x(1);-x(2:4)];
end
function P=aspoly(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
function tf=isAuto(x)
tf=(ischar(x)||(isstring(x)&&isscalar(x)))&&strcmpi(char(x),'auto');
end

function [r,t,gap,gapGuard]=svdcomparison(x0,y0,x1,y1,scaleRatio,side)
A=chiscalar(x0)+1i*chiscalar(y0);B=chiscalar(x1)+1i*chiscalar(y1);
if strcmp(side,'left'),G=B/A;else,G=A\B;end
[W,D,V]=svd(G);s=diag(D);r=scaleRatio/s(1);
gap=max(0,(s(1)-s(2))/s(1));
% A forward sensitivity estimate for the solve and 2-by-2 SVD. This is a
% numerical guard, not a validated singular-value enclosure.
gapGuard=256*eps*(1+cond(A));
kv=[real(V(1,1));imag(V(1,1));-real(V(2,1));imag(V(2,1))];
kw=[real(W(1,1));imag(W(1,1));-real(W(2,1));imag(W(2,1))];
t=-r*mulq(kv,conjq(kw));
end

function A=chiscalar(q)
a=q(1)+1i*q(2);b=q(3)+1i*q(4);A=[a,b;-conj(b),conj(a)];
end

function [q,ok]=commonorientation(C)
% A fixed quaternion times a real polynomial gives a repeated singular
% value STRUCTURALLY. Do not infer it from a small computed SVD gap.
q=[1;0;0;0];ok=false;
[~,idx]=max(max(abs(C),[],1));v=C(:,idx);[~,pivot]=max(abs(v));
if v(pivot)==0,return;end
for k=1:size(C,2)
    for j=1:4
        if j==pivot,continue;end
        if ~exactproducts(C(j,k),v(pivot),C(pivot,k),v(j)),return;end
    end
end
q=v/norm(v);ok=all(isfinite(q));
end

function yes=exactproducts(a,b,c,d)
% Error-free two-product comparison for finite binary floating-point data.
% Underflow/overflow makes the structural test inconclusive, never true.
z1=a==0||b==0;z2=c==0||d==0;
if z1||z2,yes=z1&&z2;return;end
p=a*b;q=c*d;
if ~isfinite(p)||~isfinite(q)||min(abs([p q]))<pow2(realmin,28)
    yes=false;return
end
if p~=q,yes=false;return;end
e=producterror(a,b,p);f=producterror(c,d,q);
yes=isfinite(e)&&isfinite(f)&&e==f;
end

function e=producterror(a,b,p)
split=134217729;
x=split*a;ah=x-(x-a);al=a-ah;
x=split*b;bh=x-(x-b);bl=b-bh;
e=((ah*bh-p)+ah*bl+al*bh)+al*bl;
end

function G=emptygeometry(kind,reason)
G=struct('kind',kind,'centerComponents',NaN(4,1),'basisComponents',zeros(4,0), ...
    'sphereRadius',NaN,'dimension',NaN,'selectedParameterComponents',NaN(4,1), ...
    'singularValueGapRatio',NaN,'inverseGapEstimate',NaN, ...
    'singularGapUncertaintyEstimate',NaN, ...
    'orthogonalityResidual',NaN,'boundaryResidual',NaN,'selectedMembershipResidual',NaN, ...
    'classificationReason',reason,'assurance','floating-point geometry; no exact certificate');
end

function G=gaingeometry(t,x0,y0,x1,y1,q0,q1,common,ratio,side,omega,gap,gapGuard,tol)
G=emptygeometry('undetermined','A small numerical gap does not establish exact multiplicity.');
G.selectedParameterComponents=t;G.singularValueGapRatio=gap;
G.singularGapUncertaintyEstimate=gapGuard;
if gap>0,G.inverseGapEstimate=1/gap;else,G.inverseGapEstimate=Inf;end
if omega==0
    G.kind='point';G.centerComponents=t;G.sphereRadius=0;G.dimension=0;
    G.classificationReason='At zero frequency the nonzero quaternion perturbation gives a unique critical parameter.';
elseif common
    a0=q0'*x0;b0=q0'*y0;a1=q1'*x1;b1=q1'*y1;
    den=a1*a1+b1*b1;A=(a0*a1+b0*b1)/den;
    numerator=b0*a1-a0*b1;B=numerator/den;
    if strcmp(side,'left'),left=q0;right=conjq(q1);else,left=conjq(q1);right=q0;end
    G.centerComponents=-ratio*A*mulq(left,right);
    G.sphereRadius=abs(ratio*B);
    if numerator==0&&exactproducts(b0,a1,a0,b1)
        G.kind='point';G.dimension=0;G.sphereRadius=0;
        G.classificationReason='Each polynomial has one exact coefficient orientation; the frequency quotient has zero imaginary component.';
    elseif abs(numerator)>128*eps*(abs(b0*a1)+abs(a0*b1))
        basis=zeros(4,3);
        for j=1:3,v=zeros(4,1);v(j+1)=1;basis(:,j)=mulq(mulq(left,v),right);end
        [basis,~]=qr(basis,0);
        G.kind='sphere';G.dimension=2;G.basisComponents=basis;
        G.orthogonalityResidual=norm(basis'*basis-eye(3),'fro');
        G.classificationReason='Exact coefficient collinearity supplies the repeated-value structure; the imaginary quotient is separated from zero.';
        G.inverseGapEstimate=Inf;
    else
        G.classificationReason='The repeated-value structure is known, but point versus a very small sphere is numerically unresolved.';
    end
elseif gap>max(1024*eps,gapGuard)
    G.kind='point';G.centerComponents=t;G.sphereRadius=0;G.dimension=0;
    G.classificationReason='The largest ordered-transfer singular value is numerically separated; conditioning is recorded separately.';
end
if strcmp(G.kind,'point')
    gains=G.centerComponents;G.selectedMembershipResidual=norm(t-gains)/max(realmin,norm(t));
elseif strcmp(G.kind,'sphere')
    gains=bsxfun(@plus,G.centerComponents,G.sphereRadius*[G.basisComponents,-G.basisComponents]);
    offset=t-G.centerComponents;projection=G.basisComponents'*offset;
    G.selectedMembershipResidual=max(norm(offset-G.basisComponents*projection), ...
        abs(norm(projection)-G.sphereRadius))/max(realmin,norm(t));
else
    return
end
A=chiscalar(x0)+1i*chiscalar(y0);B=chiscalar(x1)+1i*chiscalar(y1);
res=0;
for j=1:size(gains,2)
    Q=chiscalar(gains(:,j)/ratio);
    if strcmp(side,'left'),K=A+Q*B;else,K=A+B*Q;end
    s=svd(K);res=max(res,min(s)/max(realmin,norm(A,2)+norm(Q,2)*norm(B,2)));
end
G.boundaryResidual=res;
if res>100*tol||G.selectedMembershipResidual>100*tol
    G.kind='undetermined';G.dimension=NaN;
    G.classificationReason='The proposed geometry did not pass its numerical boundary or membership check.';
end
end

function Z=emptyzeros()
Z=struct('status','undetermined','rightEvaluationComponents',NaN(4,1),'leftEvaluationComponents',NaN(4,1), ...
    'rightResidual',NaN,'leftResidual',NaN,'boundaryClassResidual',NaN, ...
    'source','central quadratic remainder','assurance','floating-point', ...
    'message','No reliable pair of evaluation witnesses has been reconstructed.');
end

function Z=zerowitnesses(C0,C1,t,w,originalW,side,tol)
Z=emptyzeros();n=max(size(C0,2),size(C1,2));C0(:,end+1:n)=0;C1(:,end+1:n)=0;
C=C0;
for j=1:n
    if strcmp(side,'left'),C(:,j)=C(:,j)+mulq(t,C1(:,j));else,C(:,j)=C(:,j)+mulq(C1(:,j),t);end
end
peak=max(abs(C(:)));if peak>0,C=C/peak;end
if w==0
    Z.rightEvaluationComponents=zeros(4,1);Z.leftEvaluationComponents=zeros(4,1);
    Z.rightResidual=norm(C(:,1))/max(realmin,norm(C,'fro'));Z.leftResidual=Z.rightResidual;
    Z.boundaryClassResidual=0;
else
    [A,E]=evenodd(C,w);scale=max(abs([A;E]));
    if scale==0
        Z.status='sphericalCandidate';Z.message='The computed central quadratic remainder vanishes; exact divisibility requires an exact certificate.';return
    end
    A=A/scale;E=E/scale;den=E'*E;
    if den<=eps^2,return;end
    vr=-mulq(conjq(E),A)/den;vl=-mulq(A,conjq(E))/den;
    Z.rightEvaluationComponents=originalW*vr;Z.leftEvaluationComponents=originalW*vl;
    Z.rightResidual=evaluationresidual(C,w*vr,'right');
    Z.leftResidual=evaluationresidual(C,w*vl,'left');
    Z.boundaryClassResidual=max([abs(vr(1)),abs(vl(1)),abs(norm(vr)-1),abs(norm(vl)-1)]);
end
if max([Z.rightResidual,Z.leftResidual,Z.boundaryClassResidual])<=100*tol
    Z.status='computed';Z.message='Both evaluation-side zeros pass original-polynomial residual and boundary-class checks.';
end
end

function r=evaluationresidual(C,q,side)
y=zeros(4,1);b=0;m=norm(q);
for j=size(C,2):-1:1
    if strcmp(side,'right'),y=mulq(y,q)+C(:,j);else,y=mulq(q,y)+C(:,j);end
    b=b*m+norm(C(:,j));
end
r=norm(y)/max(realmin,b);
end
function o=parseopts(o,varargin)
    if mod(numel(varargin),2)~=0, error('qstab:freqmargin:Options','Use name/value pairs.'); end
    keys=fieldnames(o);
    for k=1:2:numel(varargin)
        if ~(ischar(varargin{k})||(isstring(varargin{k})&&isscalar(varargin{k})))
            error('qstab:freqmargin:Options','Option names must be text.');
        end
        j=find(strcmpi(char(varargin{k}),keys),1);
        if isempty(j), error('qstab:freqmargin:Options','Unknown option %s.',char(varargin{k})); end
        o.(keys{j})=varargin{k+1};
    end
end
