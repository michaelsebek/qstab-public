function [radius,info] = boxfreqmargin(P,widths,omega,varargin)
%BOXFREQMARGIN Frequency margin for one uncertain Hamilton component.
%   radius = qstab.boxfreqmargin(P,W,omega)
%   [radius,Info] = qstab.boxfreqmargin(P,W,omega,'Witness',true)
%   [radius,Info] = qstab.boxfreqmargin(P,W,omega,'Grade',g)
%
%   P is a scalar polynomial with ascending coefficients, supplied as a
%   numeric/native quaternion vector or an ordinary a,b,c,d component struct.
%   W is a finite nonnegative REAL 4-by-(g+1) array in the order real,i,j,k.
%   At most ONE row of W may contain nonzero entries. Each active coefficient
%   varies independently: abs(delta(nu,k)) <= radius*W(nu,k). Zero widths
%   fix their coefficients. An all-zero W is also accepted.
%   omega is a nonempty finite nonnegative real scalar or vector, in the
%   polynomial frequency-variable units (rad/s when s is physical time).
%
%   radius has omega's shape. Each entry minimizes the box multiplier over
%   the WHOLE imaginary quaternion sphere of radius omega(k), using the two
%   moving exclusion points of the one-component interval theorem. This is
%   a frequency result; it does not minimize over frequency, test nominal
%   stability, or enforce preservation of polynomial degree. Leading widths
%   are permitted. Use a separate degree-loss check for a global reserve.
%
%   If f_mu is component mu of P evaluated at the ordinary complex i*omega,
%   the two required values of the uncertain component are
%       zeta = +/-sqrt(-sum_{mu~=nu} f_mu^2).
%   The squares are COMPLEX SQUARES, not squared absolute values. For either
%   target the required multiplier is the maximum of its real displacement
%   divided by the total even width and its imaginary displacement divided
%   by the total odd width. A missing parity must match exactly; Tolerance
%   never projects a nonzero frozen displacement to zero.
%
%   Info.valid marks numerically resolved entries. Inf with valid=true means
%   no permitted finite perturbation gives a contact at that frequency.
%   NaN with valid=false retains unresolved/unrepresentable branch details.
%   Info.branchRadii, selectedBranch, frequencyData and witnesses retain the
%   two branches, gaps and reconstructed actual right-root witnesses.
%   Info.lower/upper are guarded floating-point estimates, not validated
%   interval arithmetic. A finite witness always supplies an actual member
%   of the original coefficient box; complex representatives are not
%   assumed to be point-evaluation zeros.
%
%   Options (name/value pairs):
%     Grade       []: infer coefficient grade; otherwise pad exact zeros.
%     Tolerance   1e-10: positive, less than one; reconstruction checks only.
%     Witness     false: true retains critical polynomials and right zeros.
%     Variable    's': this formula supports continuous time only.
%
%   Example: p=(s+1)^2+3*i*s+3*j, uncertainty only in the i*s coefficient.
%     P = struct('a',[1 2 1],'b',[0 3 0],'c',[3 0 0],'d',[0 0 0]);
%     W = zeros(4,3); W(2,2) = 1;
%     [r,I] = qstab.boxfreqmargin(P,W,1,'Witness',true);
%     r                         % 3-sqrt(5)
%     I.witnesses{1}.boundaryComponents
%
%   See also qstab.boxvalue, qstab.boxradius, qstab.boxstable, qstab.normpoly.
try
    ip=inputParser;
    addParameter(ip,'Grade',[]);
    addParameter(ip,'Tolerance',1e-10,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<1);
    addParameter(ip,'Witness',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'Variable','s',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
    parse(ip,varargin{:});o=ip.Results;
    if ~strcmpi(char(o.Variable),'s')
        error('qstab:boxfreqmargin:Variable','boxfreqmargin supports Variable=''s'' only.');
    end
    if ~(isnumeric(widths)&&isreal(widths))
        error('qstab:boxfreqmargin:Weights','W must be an explicit nonnegative real 4-by-(grade+1) array.');
    end
    if ~(isnumeric(omega)&&isreal(omega)&&isvector(omega)&&~isempty(omega)&&all(isfinite(omega(:)))&&all(omega(:)>=0))
        error('qstab:boxfreqmargin:Frequency','omega must be a nonempty finite nonnegative real scalar or vector.');
    end
    M=qstab.internal.boxmodel(P,widths,o.Grade);
    active=find(any(M.weights>0,2));
    if numel(active)>1
        error('qstab:boxfreqmargin:Components','Only one Hamilton component may be uncertain; use boxvalue for an actual quaternion point in a general box.');
    end
    if isempty(active),component=0;else,component=active;end
    radius=NaN(size(omega));nw=numel(omega);
    info=struct('schema','qstab.boxfreqmargin.1','variable','s', ...
        'scope','imaginaryQuaternionClasses','grade',M.grade,'widths',M.weights, ...
        'componentIndex',component,'omega',double(omega),'valid',false(size(omega)), ...
        'lower',zeros(size(omega)),'upper',Inf(size(omega)), ...
        'branchRadii',NaN(2,nw),'selectedBranch',zeros(size(omega)), ...
        'frequencyData',{cell(size(omega))},'witnesses',{cell(size(omega))}, ...
        'messages',{cell(size(omega))},'globalMinimumClaim',false,'validated',false, ...
        'assurance','floating-point two-point formula and reconstructed right-root checks');
    for j=1:nw
        E=localFrequency(M,component,double(omega(j)),double(o.Tolerance),o.Witness);
        radius(j)=E.radius;info.valid(j)=E.valid;info.lower(j)=E.lower;info.upper(j)=E.upper;
        info.selectedBranch(j)=E.selectedBranch;info.branchRadii(:,j)=[E.branches.radius].';
        info.witnesses{j}=E.witness;info.messages{j}=E.message;info.frequencyData{j}=E;
    end
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'),throwAsCaller(problem);end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end

function E=localFrequency(M,component,w,tol,keepWitness)
B=struct('radius',NaN,'lower',0,'upper',Inf,'valid',false,'status','undetermined', ...
    'scaledTarget',NaN,'scaledGap',[NaN;NaN],'parityRadii',[NaN;NaN], ...
    'boundaryComponents',NaN(4,1),'evaluationResidual',NaN,'witness',[],'message','');
E=struct('radius',NaN,'lower',0,'upper',Inf,'valid',false,'status','undetermined', ...
    'omega',w,'selectedBranch',0,'branches',repmat(B,1,2),'witness',[], ...
    'scaledCenter',NaN,'scaledWidths',[0;0],'scaledComponentValues',NaN(4,1), ...
    'evaluationScaleLog',NaN,'powerScaleLog',0,'message','');
n=M.grade;
if w==0,mag=[1 zeros(1,n)];
elseif w<=1,mag=w.^(0:n);
else,mag=(1/w).^(n:-1:0);E.powerScaleLog=n*log(w);end
if w>0&&any(mag==0 & (any(M.Cscaled~=0,1)|any(M.Wscaled>0,1)))
    E.message='An active monomial underflowed; no missing-parity conclusion is made.';return
end
[f,amplitude,ok]=localValue(M.Cscaled,mag);
T=bsxfun(@times,M.Wscaled,mag);widthActive=bsxfun(@and,M.Wscaled>0,mag>0);
if ~ok||any(T(:)==0 & widthActive(:))
    E.message='A nonzero value or width contribution underflowed.';return
end
if any(~isfinite([f;T(:)]))
    E.message='Evaluation overflowed.';return
end
if component==0,[~,nu]=max(abs(f));else,nu=component;end
we=sum(T(nu,1:2:end));wo=sum(T(nu,2:2:end));
E.scaledWidths=[we;wo];E.scaledCenter=f(nu);E.scaledComponentValues=f;
if M.coefficientScale>0,E.evaluationScaleLog=log(M.coefficientScale)+E.powerScaleLog;
else,E.evaluationScaleLog=-Inf;end
other=f;other(nu)=[];sc=max(abs(other));
if sc==0,zeta=0;targetError=0;
else
    z=other/sc;squares=z.^2;
    if any(other~=0 & z==0)||any(z~=0 & squares==0)
        E.message='The exclusion-point calculation loses a nonzero component by underflow.';return
    end
    h=-sum(squares);root=sqrt(h);zeta=sc*root;
    if ~isfinite(zeta)||(root~=0&&zeta==0)
        E.message='The exclusion point is outside the representable range.';return
    end
    % Near a multiple square root, use the square-root sensitivity, not a
    % user tolerance, when deciding whether a frozen parity is unattainable.
    dh=64*eps*sum(abs(squares));
    targetError=sc*(64*eps*abs(root)+sqrt(dh));
    if abs(root)>sqrt(dh),targetError=sc*(64*eps*abs(root)+dh/abs(root));end
end
err=128*eps*(n+1)*max(amplitude(nu),realmin)+targetError;
for b=1:2
    D=B;target=(3-2*b)*zeta;gap=[real(target-f(nu));imag(target-f(nu))];
    D.scaledTarget=target;D.scaledGap=gap;
    [re,se]=localRatio(gap(1),we,M,err);
    [ro,so]=localRatio(gap(2),wo,M,err);
    D.parityRadii=[re;ro];
    if strcmp(se,'infeasible')||strcmp(so,'infeasible')
        D.radius=Inf;D.lower=Inf;D.upper=Inf;D.valid=true;D.status='infeasible';
        D.message='A structurally fixed parity cannot reach this exclusion point.';
    elseif ~strcmp(se,'computed')||~strcmp(so,'computed')
        D.message='A frozen parity is within arithmetic uncertainty, or a positive radius is unrepresentable.';
    else
        D.radius=max(re,ro);
        [V,msg]=localWitness(M,nu,mag,w,gap,[we;wo],D.radius,tol);
        if isempty(V)
            D.radius=NaN;D.message=msg;
        else
            guard=max(4096*eps*(n+1),10*tol);
            D.lower=max(0,D.radius*(1-guard));D.upper=max(D.radius,V.weightedNorm)*(1+guard);
            if D.radius==0,D.upper=0;end
            D.valid=true;D.status='computed';D.boundaryComponents=V.boundaryComponents;
            D.evaluationResidual=V.evaluationResidual;
            if keepWitness,D.witness=V;end
            D.message='Two-point branch with a checked original-coefficient right-root witness.';
        end
    end
    E.branches(b)=D;
end
valid=[E.branches.valid];rr=[E.branches.radius];
if all(valid)
    [E.radius,E.selectedBranch]=min(rr);E.valid=true;
    E.lower=min([E.branches.lower]);E.upper=min([E.branches.upper]);
    if isinf(E.radius),E.status='infeasible';E.message='Neither exclusion point is attainable at this frequency.';
    else
        E.status='computed';E.witness=E.branches(E.selectedBranch).witness;
        E.message='Minimum of both exclusion-point branches at this frequency; no global frequency claim.';
    end
else
    finite=find(valid&isfinite(rr));
    if ~isempty(finite)
        [~,j]=min(rr(finite));E.selectedBranch=finite(j);E.upper=E.branches(finite(j)).upper;
        E.witness=E.branches(finite(j)).witness;E.status='bounded';
    end
    E.message='At least one exclusion-point branch is unresolved; retained finite witnesses give upper candidates only.';
end
end

function [r,status]=localRatio(gap,width,M,errorBound)
if width==0
    if gap==0,r=0;status='computed';
    elseif abs(gap)>errorBound,r=Inf;status='infeasible';
    else,r=NaN;status='undetermined';end
    return
end
if gap==0,r=0;status='computed';return;end
if M.coefficientScale==0,r=0;status='computed';return;end
r=exp(log(abs(gap))-log(width)+M.logScaleRatio);
if ~isfinite(r)||r==0,r=NaN;status='unrepresentable';else,status='computed';end
end

function [V,message]=localWitness(M,nu,mag,w,gap,widths,r,tol)
V=[];message='';D=zeros(size(M.C));
for parity=0:1
    if gap(parity+1)==0,continue;end
    [coordinate,status]=localRatio(gap(parity+1),widths(parity+1),M,0);
    if ~strcmp(status,'computed'),message='A parity perturbation is not representable.';return;end
    coordinate=sign(gap(parity+1))*coordinate;
    for k=(1+parity):2:size(M.C,2)
        if M.weights(nu,k)==0||mag(k)==0,continue;end
        D(nu,k)=coordinate*M.weights(nu,k)*(-1)^floor((k-1)/2);
        if ~isfinite(D(nu,k))||(coordinate~=0&&D(nu,k)==0)
            message='An original-coefficient perturbation overflowed or underflowed.';return
        end
    end
end
Ccrit=M.C+D;
% At zero frequency the only possible contact is q_0=0. Preserve that
% exact coefficient cancellation instead of retaining log/exp roundoff.
if w==0&&M.weights(nu,1)>0,Ccrit(nu,1)=0;end
actual=Ccrit-M.C;
if any(~isfinite([Ccrit(:);actual(:)]))
    message='Critical polynomial coefficients are not representable.';return
end
mask=M.weights>0;ratios=zeros(size(actual));ratios(mask)=abs(actual(mask))./M.weights(mask);
achieved=max(ratios(:));
if ~isfinite(achieved)||(r>0&&achieved==0)||any(actual(~mask)~=0)
    message='The reconstructed original coefficients lost the required perturbation.';return
end
if abs(achieved-r)>max(4096*eps*(M.grade+1),10*tol)*max(r,realmin)
    message='The returned coefficient budget disagrees with the two-point formula.';return
end
cs=max(abs(Ccrit(:)));
if cs==0,cc=Ccrit;else,cc=Ccrit/cs;end
if any(Ccrit(:)~=0 & cc(:)==0),message='Witness normalization erases a nonzero coefficient.';return;end
[f,amp,ok]=localValue(cc,mag);
if ~ok,message='Witness evaluation loses a nonzero coefficient contribution.';return;end
e=real(f);o=imag(f);scale=max(abs([e;o]));
if scale==0,u=[0;1;0;0];
else
    e=e/scale;o=o/scale;den=o'*o;
    if den==0,message='A nonzero even value has no odd cancellation direction.';return;end
    raw=-localMul([o(1);-o(2:4)],e)/den;
    guard=max(4096*eps*(M.grade+1),10*tol);
    if ~all(isfinite(raw))||abs(raw(1))>guard||abs(norm(raw)-1)>guard||norm(raw(2:4))==0
        message='The reconstructed contact does not give a unit pure right-root direction.';return
    end
    u=[0;raw(2:4)/norm(raw(2:4))];
end
res=real(f)+localMul(imag(f),u);
relativeResidual=norm(res)/max([sum(amp),realmin]);
point=w*u;
if any(~isfinite(point))||any(point==0 & u~=0 & w>0)
    message='The actual quaternion boundary point is not representable.';return
end
if relativeResidual>max(4096*eps*(M.grade+1),10*tol)
    message='The returned original polynomial failed its right-evaluation residual check.';return
end
V=struct('kind','finiteZero','valid',true,'evaluationSide','right', ...
    'boundaryComponents',point,'classRepresentative',complex(0,w), ...
    'coefficientComponents',Ccrit,'perturbationComponents',actual, ...
    'polynomial',localParts(Ccrit),'weightedNorm',achieved, ...
    'formulaRadius',r,'evaluationResidual',relativeResidual, ...
    'fixedComponentsPreserved',true,'weights',M.weights,'grade',M.grade,'validated',false);
end

function [f,amplitude,ok]=localValue(C,mag)
T=bsxfun(@times,C,mag);nz=bsxfun(@and,C~=0,mag>0);
ok=~any(T(:)==0 & nz(:));
signPower=(-1).^floor((0:size(C,2)-1)/2);
T=bsxfun(@times,T,signPower);
f=complex(sum(T(:,1:2:end),2),sum(T(:,2:2:end),2));
amplitude=sum(abs(T),2);
end
function z=localMul(x,y)
z=[x(1)*y(1)-x(2:4)'*y(2:4);x(1)*y(2:4)+y(1)*x(2:4)+cross(x(2:4),y(2:4))];
end
function P=localParts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
