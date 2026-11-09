function R = stabilityradiusReport(P0, P1, varargin)
%STABILITYRADIUSREPORT Internal implementation returning the full report.

    o=struct('Side','left','AbsTol',1e-6,'RelTol',1e-5,...
        'Tolerance',1e-10,'MaxCells',4096,'MaxTime',30,'Variable','s', ...
        'Grade',[],'MaxDepth',40,'StoreCells',true);
    o=parseopts(o,varargin{:}); side=lower(char(o.Side));
    if ~any(strcmp(side,{'left','right'})), error('qstab:stabilityradius:Side','Invalid Side.'); end
    if any(strcmpi(o.Variable,{'z','zinv'}))
        R=qstab.internal.discreteradius(P0,P1,o); return
    elseif ~strcmpi(o.Variable,'s')
        error('qstab:stabilityradius:Variable','Variable must be s, z or zinv.');
    end
    if ~isempty(o.Grade)
        error('qstab:stabilityradius:Grade','Grade is currently a discrete-radius option only.');
    end
    positive={'AbsTol','RelTol','Tolerance'};
    for j=1:numel(positive)
        v=o.(positive{j});
        if ~(isnumeric(v)&&isreal(v)&&isscalar(v)&&isfinite(v)&&v>0)
            error('qstab:stabilityradius:Options','Tolerances must be positive finite scalars.');
        end
    end
    if ~(isnumeric(o.MaxCells)&&isreal(o.MaxCells)&&isscalar(o.MaxCells)&&...
            isfinite(o.MaxCells)&&o.MaxCells>=1&&o.MaxCells==floor(o.MaxCells)) || ...
            ~(isnumeric(o.MaxTime)&&isreal(o.MaxTime)&&isscalar(o.MaxTime)&&o.MaxTime>0)
        error('qstab:stabilityradius:Budget','MaxCells must be a positive integer and MaxTime positive.');
    end
    C0=qstab.polycomponents(P0); C1=qstab.polycomponents(P1);
    n0=find(any(C0~=0,1),1,'last'); n1=find(any(C1~=0,1),1,'last');
    R=struct('status','assumptionFailed','method','global-frequency-branch-and-bound',...
        'assurance','floating-point-bounds','validated',false,'side',side,...
        'variable','s','lower',NaN,'upper',NaN,'estimate',NaN,'gap',NaN,...
        'frequency',NaN,'criticalParameterComponents',NaN(4,1),...
        'boundaryComponents',NaN(4,1),'message','','diagnostics',struct());
    if isempty(n0)||n0<2||(~isempty(n1)&&n1>=n0)
        R.message='Require degree(P0)>=1 and degree(P1)<degree(P0).'; return
    end
    nominal=qstab.internal.isstableReport(P0,'Variable','s','Tolerance',o.Tolerance); R.nominal=nominal;
    if ~strcmp(nominal.status,'stable')
        if strcmp(nominal.status,'undetermined'), R.status='undetermined'; end
        R.message='P0 must be numerically verified Hurwitz.'; return
    end
    if isempty(n1)
        R.status='unbounded'; R.lower=Inf; R.upper=Inf; R.estimate=Inf; R.gap=0;
        R.message='P1 is identically zero; the family is constant.'; return
    end
    C0=C0(:,1:n0); C1=C1(:,1:n1); degree=n0-1;
    a=hypot(hypot(C0(1,:),C0(2,:)),hypot(C0(3,:),C0(4,:)));
    b=zeros(1,degree);
    b(1:n1)=hypot(hypot(C1(1,:),C1(2,:)),hypot(C1(3,:),C1(4,:)));
    tstart=tic; count=0; upperBound=Inf; bestW=NaN; bestI=[];
    initial=unique([0,1,0:degree]);
    [initialR,initialInfo]=qstab.freqmargin(P0,P1,initial,'Side',side,...
        'CheckNominal',false,'Tolerance',o.Tolerance);
    count=count+numel(initial);
    for k=1:numel(initial)
        if initialInfo.valid(k)&&isfinite(initialR(k))&&initialR(k)<upperBound
            upperBound=initialR(k); bestW=initial(k); bestI=pickinfo(initialInfo,k);
        end
    end
    if ~isfinite(upperBound)
        R.status='undetermined'; R.lower=0;
        R.message='No reliable finite initial frequency witness was obtained.'; return
    end
    % For w >= cutoff, |a_n|*w^n - sum_{k<n}(|a_k|+U|b_k|)*w^k > 0.
    % Hence all boundary parameters there have modulus greater than U.
    cutoff=1+max((a(1:degree)+upperBound*b)/a(end));
    if ~isfinite(cutoff)
        R.status='bounded'; R.lower=0; R.upper=upperBound; R.estimate=upperBound; R.gap=upperBound;
        R.frequency=bestW; R.criticalParameterComponents=bestI.parameter;
        R.boundaryComponents=bestI.boundary;
        R.message='Frequency-tail bound overflowed; only the witness upper bound is available.'; return
    end
    tailLower=upperBound;
    % Cells store [left,right,lower,midpoint objective]. No gridded
    % acceptance: the minimum lower bound covers the whole frequency range.
    cellData=zeros(max(2,o.MaxCells+1),4);
    [cellData(1,:),rv,iv]=makecell(0,cutoff);
    ncell=1; updatebest(rv,iv,cutoff/2);
    shortcut=bernsteinlower(C0,C1,upperBound,tailLower,cutoff,side,o,tstart);
    if shortcut.accepted
        cellData(1,3)=max(cellData(1,3),shortcut.radius);
    end
    status='bounded'; reason='Cell or time budget reached.';
    while true
        lowerBound=min(tailLower,min(cellData(1:ncell,3)));
        if lowerBound>upperBound
            % Roundoff cannot be silently hidden as a mathematically valid bracket.
            if lowerBound-upperBound>100*eps*max(1,upperBound)
                status='undetermined'; reason='Computed lower and upper bounds are inconsistent.'; break
            end
            lowerBound=max(0,upperBound-100*eps*max(1,upperBound));
        end
        if upperBound-lowerBound<=o.AbsTol+o.RelTol*upperBound
            status='converged'; reason='Global numerical bound gap meets the requested tolerance.';
            if shortcut.accepted
                reason='A guarded Bernstein norm inequality supplies the global numerical lower bound.';
            end
            break
        end
        if ncell>=o.MaxCells||toc(tstart)>=o.MaxTime, break; end
        [~,idx]=min(cellData(1:ncell,3)); cellNow=cellData(idx,:);
        l=cellNow(1); r=cellNow(2); parentMidpoint=l+(r-l)/2;
        if parentMidpoint==l||parentMidpoint==r
            reason='Frequency subdivision reached floating-point resolution.'; break
        end
        [left,rv,iv]=makecell(l,parentMidpoint); updatebest(rv,iv,l+(parentMidpoint-l)/2);
        [right,rv,iv]=makecell(parentMidpoint,r); updatebest(rv,iv,parentMidpoint+(r-parentMidpoint)/2);
        cellData(idx,:)=left; ncell=ncell+1; cellData(ncell,:)=right;
    end
    R.status=status; R.lower=lowerBound; R.upper=upperBound; R.estimate=upperBound; R.gap=upperBound-lowerBound;
    R.frequency=bestW; R.criticalParameterComponents=bestI.parameter;
    R.boundaryComponents=bestI.boundary; R.message=reason;
    R.diagnostics=struct('frequencyCutoff',cutoff,'tailLowerBound',tailLower,...
        'cells',cellData(1:ncell,:),'cellColumns',{{'left','right','lower','midpointValue'}},...
        'numberOfCells',ncell,'frequencyEvaluations',count,'elapsedSeconds',toc(tstart),...
        'witnessEvaluationResidual',bestI.residual,'AbsTol',o.AbsTol,'RelTol',o.RelTol);
    R.diagnostics.bernsteinShortcut=shortcut;
    R.diagnostics.bernsteinCells=shortcut.processedCells;
    R.diagnostics.lowerBoundMethod='frequency-cell-lipschitz';
    if shortcut.accepted
        R.diagnostics.lowerBoundMethod='bernstein-norm-inequality';
    end

    function [row,rv,iv]=makecell(l,r)
        m=l+(r-l)/2; h=(r-l)/2;
        [rv,iv]=qstab.freqmargin(P0,P1,m,'Side',side,'CheckNominal',false,'Tolerance',o.Tolerance);
        count=count+1; lb=0;
        if iv.valid
            mn=iv.minNominalModulus*(1-100*o.Tolerance);
            da=h*derivbound(a,r); db=h*derivbound(b,r);
            if isfinite(da)&&isfinite(db)&&mn>da
                den=mn/iv.centerLowerEstimate+db;
                if den==0, lb=Inf; else, lb=(mn-da)/den; end
                if isfinite(lb), lb=max(0,lb-100*eps*max(1,lb)); end
            end
        end
        row=[l,r,lb,rv];
    end
    function updatebest(rv,iv,w)
        if iv.valid&&isfinite(rv)&&rv<upperBound
            upperBound=rv; bestW=w; bestI=pickinfo(iv,1);
        end
    end
end

function S=bernsteinlower(C0,C1,U,tail,W,side,o,timer)
% With x=w^2, N-r^2*D >= h(x)-2*sqrt(x)*|d(x)|. The two polynomial
% inequalities h>=0 and g=h^2-4*x*|d|^2>=0 therefore imply rho(w)>=r.
% Per-coefficient error estimates below are conservative numerical guards,
% not outward-rounded interval arithmetic. No negative coefficient is clipped.
    r=max(0,U-.9*(o.AbsTol+o.RelTol*U));
    S=struct('method','bernstein-norm-inequality','attempted',true,'accepted',false,...
        'radius',r,'processedCells',0,'acceptedCells',0,'maximumDepth',0,...
        'minimumGuardedCoefficient',Inf,'validated',false,'reason','Positivity not established.');
    if ~isfinite(r)||r>tail, S.reason='Target is incompatible with the tail bound.'; return; end
    if r==0, S.accepted=true; S.reason='Nonnegativity of the radius gives the zero lower bound.'; return; end
    X=W*W;
    s0=max(abs(C0(:))); s1=max(abs(C1(:)));
    ratio=exp(log(r)+log(s1)-log(s0));
    if ~isfinite(X)||~isfinite(ratio)||ratio==0
        S.reason='Squared cutoff or relative coefficient scale is unrepresentable.'; return
    end
    common=max(1,ratio);
    C0=(C0/s0)/common; C1=(C1/s1)*(ratio/common);
    [n0,d0,a0,b0]=normparts(C0,side); [n1,d1,a1,b1]=normparts(C1,side);
    h=padd(n0,-n1); ha=padd(a0,a1);
    count=max(size(d0,2),size(d1,2));
    d0(:,end+1:count)=0; d1(:,end+1:count)=0;
    b0(:,end+1:count)=0; b1(:,end+1:count)=0;
    [ds,~]=psquares(d0-d1); [da,~]=psquares(b0+b1);
    g=padd(conv(h,h),-4*[0,ds]); ga=padd(conv(ha,ha),4*[0,da]);
    gamma=256*eps*(numel(g)+1)^2;
    [bh,eh]=powerbernstein(h,gamma*ha,X);
    [bg,eg]=powerbernstein(g,gamma*ga,X);
    if any(~isfinite([bh,eh,bg,eg]))
        S.reason='Polynomial or Bernstein conversion overflowed.'; return
    end
    stack={{bh,eh,bg,eg,0}}; limit=min(128,o.MaxCells);
    while ~isempty(stack)
        if S.processedCells>=limit||toc(timer)>=o.MaxTime
            S.reason='Bernstein shortcut budget reached; frequency subdivision remains available.'; return
        end
        cellNow=stack{end}; stack(end)=[];
        S.processedCells=S.processedCells+1;
        S.maximumDepth=max(S.maximumDepth,cellNow{5});
        guarded=[cellNow{1}-cellNow{2},cellNow{3}-cellNow{4}];
        if all(guarded>0)
            S.acceptedCells=S.acceptedCells+1;
            S.minimumGuardedCoefficient=min(S.minimumGuardedCoefficient,min(guarded));
            continue
        end
        if cellNow{5}>=30 || max(cellNow{1}+cellNow{2})<0 || max(cellNow{3}+cellNow{4})<0
            S.reason='A polynomial sign or arithmetic guard prevents this sufficient test.'; return
        end
        [hl,el,hr,er]=splitbernstein(cellNow{1},cellNow{2});
        [gl,fl,gr,fr]=splitbernstein(cellNow{3},cellNow{4});
        stack{end+1}={hr,er,gr,fr,cellNow{5}+1}; %#ok<AGROW>
        stack{end+1}={hl,el,gl,fl,cellNow{5}+1}; %#ok<AGROW>
    end
    S.accepted=true; S.reason='All covering Bernstein cells satisfy both guarded polynomial inequalities.';
end

function [n,d,na,da]=normparts(C,side)
    e=bsxfun(@times,C(:,1:2:end),(-1).^(0:size(C(:,1:2:end),2)-1));
    odd=C(:,2:2:end);
    if isempty(odd), odd=zeros(4,1); end
    odd=bsxfun(@times,odd,(-1).^(0:size(odd,2)-1));
    [ne,ae]=psquares(e); [no,ao]=psquares(odd);
    n=padd(ne,[0,no]); na=padd(ae,[0,ao]);
    d=zeros(3,size(e,2)+size(odd,2)-1); da=d;
    crossSign=-1; if strcmp(side,'right'), crossSign=1; end
    for k=1:3
        j=mod(k,3)+1; ell=mod(k+1,3)+1;
        d(k,:)=conv(odd(1,:),e(k+1,:))-conv(odd(k+1,:),e(1,:))+...
            crossSign*(conv(odd(j+1,:),e(ell+1,:))-conv(odd(ell+1,:),e(j+1,:)));
        da(k,:)=conv(abs(odd(1,:)),abs(e(k+1,:)))+conv(abs(odd(k+1,:)),abs(e(1,:)))+...
            conv(abs(odd(j+1,:)),abs(e(ell+1,:)))+conv(abs(odd(ell+1,:)),abs(e(j+1,:)));
    end
end
function [p,a]=psquares(C)
    p=zeros(1,2*size(C,2)-1); a=p;
    for k=1:size(C,1)
        p=p+conv(C(k,:),C(k,:)); a=a+conv(abs(C(k,:)),abs(C(k,:)));
    end
end
function p=padd(a,b)
    n=max(numel(a),numel(b)); a(end+1:n)=0; b(end+1:n)=0; p=a+b;
end
function [b,e]=powerbernstein(p,pe,X)
    n=numel(p)-1; powers=[1,cumprod(X*ones(1,n))];
    p=p.*powers; pe=pe.*powers; b=zeros(1,n+1); e=b;
    for k=0:n
        ratio=1; absoluteSum=0;
        for j=0:k
            if j>0, ratio=ratio*(k-j+1)/(n-j+1); end
            b(k+1)=b(k+1)+ratio*p(j+1);
            e(k+1)=e(k+1)+ratio*pe(j+1);
            absoluteSum=absoluteSum+abs(ratio*p(j+1));
        end
        e(k+1)=e(k+1)+16*eps*(n+1)*absoluteSum;
    end
end
function [l,le,r,re]=splitbernstein(b,e)
    n=numel(b); l=zeros(size(b)); r=l; le=l; re=l;
    l(1)=b(1); r(n)=b(end); le(1)=e(1); re(n)=e(end);
    for k=1:n-1
        e=.5*e(1:end-1)+.5*e(2:end)+4*eps*(.5*abs(b(1:end-1))+.5*abs(b(2:end)));
        b=.5*b(1:end-1)+.5*b(2:end);
        l(k+1)=b(1); r(n-k)=b(end); le(k+1)=e(1); re(n-k)=e(end);
    end
end

function v=derivbound(coeff,w)
    v=0;
    for k=numel(coeff)-1:-1:1, v=v*w+k*coeff(k+1); end
end
function I=pickinfo(info,k)
    I=struct('parameter',info.parameterComponents(:,k),...
        'boundary',info.boundaryComponents(:,k),'residual',info.evaluationResidual(k));
end
function o=parseopts(o,varargin)
    if mod(numel(varargin),2)~=0, error('qstab:stabilityradius:Options','Use name/value pairs.'); end
    keys=fieldnames(o);
    for k=1:2:numel(varargin)
        if ~(ischar(varargin{k})||(isstring(varargin{k})&&isscalar(varargin{k})))
            error('qstab:stabilityradius:Options','Option names must be text.');
        end
        j=find(strcmpi(char(varargin{k}),keys),1);
        if isempty(j), error('qstab:stabilityradius:Options','Unknown option %s.',char(varargin{k})); end
        o.(keys{j})=varargin{k+1};
    end
end
