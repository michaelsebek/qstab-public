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
%     CheckNominal true; false skips the nominal stability check when its
%                  result is already known.
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
    o=struct('Side','left','Tolerance',1e-10,'CheckNominal',true,'Variable','s','Grade',[]);
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
    if o.CheckNominal
        nominal=qstab.internal.isstableReport(P0,'Variable','s','Tolerance',o.Tolerance);
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
        'nominal',nominal,'message',{cell(size(omega))});
    for k=1:nw
        w=double(omega(k));
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
        info.boundaryComponents(:,k)=w*vq;
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
    end
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
