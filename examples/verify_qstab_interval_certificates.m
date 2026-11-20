function report = verify_qstab_interval_certificates(varargin)
%VERIFY_QSTAB_INTERVAL_CERTIFICATES Recompute the article's exact certificates.
%   REPORT = VERIFY_QSTAB_INTERVAL_CERTIFICATES recomputes all 56 quartic
%   Bernstein inequalities and the complete boundaries of all four explicit
%   quintics and their central Routh counts, using Symbolic Math Toolbox
%   and exact rational arithmetic.
%
%   REPORT = VERIFY_QSTAB_INTERVAL_CERTIFICATES('Part','quartic') recomputes
%   only the quartic inequalities. 'Part','quintic' checks only the quintics.
%   REPORT = VERIFY_QSTAB_INTERVAL_CERTIFICATES(...,'Verbose',false) suppresses
%   progress printing. The default Part is 'all', and Verbose is true.
%
%   REPORT.ok is true only after all requested claims were recomputed in this
%   invocation. REPORT.computedThisRun explicitly distinguishes this result
%   from supplied external certificates. This function does not read stored
%   pass flags as evidence. It creates no files and changes no path or folder.
%
%   The source quaternion polynomial components are constructed from exact
%   rational constants. For each quintic boundary face, the real norm and its
%   first nine Hurwitz determinants are recomputed symbolically. Their whole
%   two-dimensional domains are covered by retained dyadic Bernstein cells;
%   every coefficient is compared exactly, and complete coverage is checked.
%   The positive constant coefficient times determinant nine is determinant
%   ten. Thus this is a continuous-face proof, not a parameter-grid test.
%
%   Requirements: MATLAB and Symbolic Math Toolbox. This optional calculation
%   is independent of the mandatory base-MATLAB study. Symbolic determinant
%   computation can take substantially longer than the numerical study.
%   Its first native MATLAB validation is pending for this delivery.
%
%   Independent Python replay (requires SymPy) is also supplied in
%   reproducibility/qb02. Neither replay implements the article's general
%   exact fixed-box or global algebraic reserve algorithm; these are specific
%   reproducible certificates for its explicitly given examples.
%
%   Example:
%     proof = verify_qstab_interval_certificates('Part','all');
%     assert(proof.ok);
%     proof.coverage
%
%   See also run_qstab_intervals, qstab_interval_cases.

ip=inputParser;
ip.FunctionName=mfilename;
addParameter(ip,'Part','all',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'Verbose',true,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});
part=lower(char(ip.Results.Part));
if ~ismember(part,{'all','quartic','quintic'})
    error('qstab:exact:Part','Part must be all, quartic, or quintic.');
end
if isempty(which('sym'))
    error('qstab:exact:SymbolicRequired', ...
        ['This optional exact replay requires Symbolic Math Toolbox. ', ...
         'The ordinary interval study does not require it.']);
end
verbose=ip.Results.Verbose;
ticId=tic;
packageRoot=fileparts(fileparts(mfilename('fullpath')));
certDir=fullfile(packageRoot,'reproducibility','qb02');
report=struct('ok',false,'status','RUNNING','computedThisRun',true, ...
    'arithmetic','exact rational symbolic','dependency','Symbolic Math Toolbox', ...
    'part',part,'quarticInequalities',0,'quinticFaces',0, ...
    'quinticPolynomialChecks',0,'proofLeaves',0,'quinticCentres',0, ...
    'centralRouth',struct('family',{},'negativeIndexZeroBased',{},'negativeEntry',{},'signs',{},'rightHalfPlaneNormRoots',{}), ...
    'elapsedSeconds',NaN, ...
    'provenance','Recomputed from original rational coefficients; stored pass flags are not evidence.');
s=sym('s'); v=sym('v'); delta=sym('delta'); x=sym('x'); y=sym('y');
if ismember(part,{'all','quartic'})
    data=jsondecode(fileread(fullfile(certDir,'exact_example_certificate.json')));
    require(strcmp(data.schema,'qstab.interval-exact-examples.1'),'Wrong quartic certificate schema.');
    claims=quarticClaims(s,v,delta);
    require(numel(data.bernstein)==numel(claims),'The quartic certificate is incomplete.');
    for k=1:numel(claims)
        c=claims(k); saved=data.bernstein(k);
        require(strcmp(saved.name,c.name),'Quartic claim order differs.');
        require(equalQ(rational(item(saved.interval,1)),c.lo) && ...
            equalQ(rational(item(saved.interval,2)),c.hi) && ...
            equalQ(rational(saved.shift_max),c.deltaMax),'Quartic proof rectangle differs.');
        require(equalQ(rational(saved.claimed_strict_lower_bound),c.target),'Quartic bound differs.');
        P=expand(subs(c.polynomial,[v delta],[c.lo+(c.hi-c.lo)*x,c.deltaMax*y]));
        B=bernsteinUnit(P,x,y,double(saved.bidegree(1)),double(saved.bidegree(2)));
        require(all(isAlways(B(:)>c.target)),'A quartic Bernstein bound is not strictly positive.');
        verifyMinimum(B,rational(saved.minimum_exact));
        report.quarticInequalities=report.quarticInequalities+1;
        if verbose, fprintf('PASS [exact quartic] %s\n',c.name); end
    end
end
if ismember(part,{'all','quintic'})
    data=jsondecode(fileread(fullfile(certDir,'quintic_boundary_certificate.json')));
    require(strcmp(data.schema,'qstab.interval-quintic-boundaries.1'),'Wrong quintic certificate schema.');
    require(numel(data.faces)==24,'The complete four-box boundary requires 24 faces.');
    names={'OOO','EOO','EEE','EEO'}; index=0;
    for f=1:4
        centre=quinticComponents(names{f},s,sym(zeros(1,3)));
        centreNorm=expand(sum(centre.^2));
        column=exactRouth(centreNorm,s);
        if f<=2, neg=10; else, neg=8; end
        if ismember(f,[1 3])
            lower=-rational('1/500000'); upper=-rational('1/1000000');
        else
            lower=-rational('7/100000000'); upper=-rational('3/50000000');
        end
        other=setdiff(1:11,neg);
        require(all(isAlways(column(other)>0)) && isAlways(column(neg)>lower) && ...
            isAlways(column(neg)<upper),'The exact quintic central Routh bounds failed.');
        signs=ones(1,11); signs(neg)=-1;
        changes=sum(signs(1:end-1).*signs(2:end)<0);
        require(changes==2,'The central Routh column does not have two sign changes.');
        report.centralRouth(end+1)=struct('family',names{f}, ...
            'negativeIndexZeroBased',neg-1,'negativeEntry',char(column(neg)), ...
            'signs',signs,'rightHalfPlaneNormRoots',changes);
        report.quinticCentres=report.quinticCentres+1;
        if verbose, fprintf('PASS [exact centre] %s, two RHP norm roots\n',names{f}); end
        for fixed=1:3
            for sg=[-1 1]
                index=index+1; face=data.faces(index);
                require(strcmp(face.family,names{f}) && face.fixed_coordinate==fixed-1 && face.sign==sg, ...
                    'A retained boundary face differs from the original box.');
                require(numel(face.checks)==11,'A face is missing Hurwitz polynomial checks.');
                t=sym(zeros(1,3)); free=setdiff(1:3,fixed);
                t(fixed)=sym(sg); t(free)=[x y];
                C=quinticComponents(names{f},s,t);
                N=expand(sum(C.^2));
                a=coeffs(N,s,'All'); a=reshape(a,1,numel(a));
                require(numel(a)==11,'The quintic norm must have degree ten.');
                require(equalQ(expand(sum(a.*s.^(10:-1:0))-N),sym(0)),'Norm coefficient order differs.');
                for q=1:11
                    if q==1
                        P=a(1); label='leading';
                    elseif q==2
                        P=a(end); label='constant';
                    else
                        n=q-2; H=sym(zeros(n));
                        for ii=0:n-1
                            for jj=0:n-1
                                ind=2*jj-ii+1;
                                if ind>=0 && ind<=10, H(ii+1,jj+1)=a(ind+1); end
                            end
                        end
                        P=expand(det(H)); label=sprintf('Delta%d',n);
                    end
                    saved=face.checks(q);
                    require(strcmp(saved.name,label),'A Hurwitz determinant is missing or reordered.');
                    % First transform the actual whole face [-1,1]^2 to [0,1]^2.
                    P=expand(subs(P,[x y],[2*x-1,2*y-1]));
                    B=bernsteinUnit(P,x,y,double(saved.bidegree(1)),double(saved.bidegree(2)));
                    leaves=saved.leaves;
                    paths=arrayfun(@(z)char(z.path),leaves,'UniformOutput',false);
                    lower=arrayfun(@(z)char(z.lower),leaves,'UniformOutput',false);
                    replayLeaves(B,paths,lower);
                    report.quinticPolynomialChecks=report.quinticPolynomialChecks+1;
                    report.proofLeaves=report.proofLeaves+numel(leaves);
                end
                report.quinticFaces=report.quinticFaces+1;
                if verbose
                    fprintf('PASS [exact quintic] %s, t%d=%+d/1000; complete face\n',names{f},fixed,sg);
                end
            end
        end
    end
end
report.ok=true; report.status='EXACT_CERTIFICATES_RECOMPUTED';
report.elapsedSeconds=toc(ticId);
report.coverage=struct('quarticInequalities',report.quarticInequalities, ...
    'completeQuinticFaces',report.quinticFaces, ...
    'quinticPolynomialChecks',report.quinticPolynomialChecks, ...
    'dyadicProofLeaves',report.proofLeaves,'quinticCentres',report.quinticCentres);
if verbose
    fprintf('Exact interval certificates: %s | %.3f s\n',report.status,report.elapsedSeconds);
end
end

function C=baseComponents(kind,s)
if strcmp(kind,'O')
    C=[s^4+rational('1/10')*s^3+rational('1/4')*s^2-rational('1/10')*s-rational('1/8'), ...
       rational('5/4')*s^3+rational('1/5')*s^2+rational('5/8')*s, ...
       rational('3/40')*(s^2+1), sym(0)];
else
    C=[s^4+rational('3/2')*s^3-rational('1/4')*s^2-rational('1/2')*s-1, ...
       rational('17/8')*s^3+2*s^2+rational('15/8')*s, ...
       -rational('3/2')*(s^3+s), sym(0)];
end
end

function C=quinticComponents(name,s,t)
r=rational('1/1000'); ep=rational('1/1000000000000');
if ismember(name,{'OOO','EEE'})
    C=subs(baseComponents('O',s),s,s-rational('1/4000000'));
    C(2:4)=C(2:4)+r*t*s;
else
    C=subs(baseComponents('E',s),s,s-rational('1/125000000'));
    C(2)=C(2)+r*t(1); C(3:4)=C(3:4)+r*t(2:3)*s;
end
if ismember(name,{'OOO','EOO'})
    C(1)=C(1)+ep*s^5;
else
    C=s*C; C(1)=C(1)-ep;
end
C=expand(C);
end

function column=exactRouth(N,s)
a=coeffs(N,s,'All'); a=reshape(a,1,numel(a)); n=numel(a)-1; m=ceil((n+1)/2);
T=sym(zeros(n+1,m));
T(1,1:numel(a(1:2:end)))=a(1:2:end);
T(2,1:numel(a(2:2:end)))=a(2:2:end);
for row=3:n+1
    require(isAlways(T(row-1,1)~=0),'An exact Routh pivot is zero.');
    for k=1:m-1
        T(row,k)=(T(row-1,1)*T(row-2,k+1)-T(row-2,1)*T(row-1,k+1))/T(row-1,1);
    end
end
column=T(:,1);
end

function [E,O]=evenOdd(C,s,v)
E=sym(zeros(1,4)); O=E;
for k=0:4
    c=subs(diff(C,s,k),s,0)/factorial(sym(k));
    if mod(k,2)==0
        E=E+(-1)^(k/2)*c*v^(k/2);
    else
        O=O+(-1)^((k-1)/2)*c*v^((k-1)/2);
    end
end
E=expand(E); O=expand(O);
end

function C=quarticClaims(s,v,delta)
C=struct('name',{},'polynomial',{},'lo',{},'hi',{},'deltaMax',{},'target',{});
r=rational('1/1000'); deltaO=rational('1/4000000'); deltaE=rational('1/125000000');
[E,O]=evenOdd(subs(baseComponents('O',s),s,s-delta),s,v);
a=E(1); b=O(1); c1=E(2); c2=E(3); y1=O(2); y2=O(3);
D=c1^2+c2^2; P=v*b^2-D; n=[-a*b*c1-D*y1,-a*b*c2-D*y2,sym(0)];
C=addClaim(C,'OOO D-3/100',D-rational('3/100'),'9/10','11/10',deltaO,'0');
C=addClaim(C,'OOO 41/1000-D',rational('41/1000')-D,'997/1000','1003/1000',deltaO,'0');
C=addClaim(C,'OOO a2-4D',a^2-4*D,'9/10','11/10',deltaO,'0');
C=addClaim(C,'OOO P lower',P,'9/10','997/1000',deltaO,'1/1000000000');
C=addClaim(C,'OOO P upper',P,'1003/1000','11/10',deltaO,'1/1000000000');
js=[1 1 2 2 3]; signs=[1 -1 1 -1 1];
targets={'1/5000000000','1/5000000000','3/100000000000','7/1000000000000','1/1250000000'};
for k=1:numel(js)
    j=js(k); sg=signs(k); H=v*(sg*r*D-n(j))^2+D*(D+a^2)*P;
    C=addClaim(C,sprintf('OOO H%d,%+d',j,sg),H,'997/1000','1003/1000',deltaO,targets{k});
end
[E,O]=evenOdd(subs(baseComponents('E',s),s,s-delta),s,v);
a=E(1); b=O(1); u=E(2); dd=O(2); c=E(3); y=O(3);
cc=v*dd^2-c^2; P=cc-a^2; n2=a*b*c-cc*y;
B=a*b+dd*u+c*y; Ab=a^2+u^2+c^2-v*(b^2+dd^2+y^2);
C=addClaim(C,'EOO C-1/20',cc-rational('1/20'),'99/100','101/100',deltaE,'0');
C=addClaim(C,'EOO 3/40-C',rational('3/40')-cc,'99/100','101/100',deltaE,'0');
C=addClaim(C,'EOO 3/40-d2',rational('3/40')-dd^2,'99/100','101/100',deltaE,'0');
C=addClaim(C,'EOO P lower',P,'99/100','4997/5000',deltaE,'3/2000000000');
C=addClaim(C,'EOO P upper',P,'5003/5000','101/100',deltaE,'3/2000000000');
for sg=[-1 1]
    for tau=[-1 1]
        C=addClaim(C,sprintf('EOO face B %d %d',sg,tau),-sg*(B+dd*sg*r+c*tau*r), ...
            '4997/5000','5003/5000',deltaE,'17/100000');
    end
    H=cc*(sg*r*cc-n2)^2+dd^2*cc*(cc+v*b^2)*P;
    C=addClaim(C,sprintf('EOO H2,%+d',sg),H,'4997/5000','5003/5000',deltaE,'7/10000000000000');
end
H3=v*r^2*cc+(cc+v*b^2)*P;
C=addClaim(C,'EOO H3',H3,'4997/5000','5003/5000',deltaE,'1/25000000');
vi={'3/8','2/5','17/20','7/8','99/100','101/100'}; signs=[1 -1 -1 1 1 -1];
for k=1:numel(vi)
    for sg=[-1 1]
        for tau=[-1 1]
            t1=sg*r; t2=tau*r;
            name=sprintf('EOO bracket B %s %d/1000 %d/1000',vi{k},sg,tau);
            C=addClaim(C,name,subs(signs(k)*(B+dd*t1+c*t2),v,rational(vi{k})), ...
                '0','1',deltaE,'9/10000');
        end
    end
end
lo={'3/8','17/20'}; hi={'2/5','7/8'};
a0=[-rational('1/25'),rational('3/20')]; umax=[sym(1),sym(2)];
vymax=[rational('2/5'),rational('3/10')];
for k=1:2
    if k==1, p=a0(k)-Ab; else, p=Ab-a0(k); end
    C=addClaim(C,['EOO Ab ' lo{k}],p,lo{k},hi{k},deltaE,'0');
    for sg=[-1 1]
        C=addClaim(C,sprintf('EOO u %s %d',lo{k},sg),umax(k)-sg*u,lo{k},hi{k},deltaE,'0');
        C=addClaim(C,sprintf('EOO vy %s %d',lo{k},sg),vymax(k)-sg*v*y,lo{k},hi{k},deltaE,'0');
    end
end
require(numel(C)==56,'Internal quartic claim catalogue is incomplete.');
end

function C=addClaim(C,name,P,lo,hi,dm,target)
C(end+1)=struct('name',name,'polynomial',expand(P),'lo',rational(lo), ...
    'hi',rational(hi),'deltaMax',dm,'target',rational(target));
end

function B=bernsteinUnit(P,x,y,m,n)
% All coefficients are retained, including internal zeros. Read each power
% row at once rather than dispatching a symbolic operation for every term.
cx=fliplr(coeffs(expand(P),x,'All')); A=sym(zeros(m+1,n+1));
require(numel(cx)<=m+1,'Retained x degree is too small.');
for i=1:numel(cx)
    cy=fliplr(coeffs(cx(i),y,'All'));
    require(numel(cy)<=n+1,'Retained y degree is too small.');
    require(isempty(symvar(cy)),'Bernstein conversion left a symbolic coefficient.');
    A(i,1:numel(cy))=cy;
end
B=bernsteinTransform(m)*A*bernsteinTransform(n).';
end

function T=bernsteinTransform(m)
% Integer numerator/denominator arrays become symbolic BEFORE division.
% Cache only the universal basis change, never a claim or proof outcome.
persistent orders transforms
if isempty(orders),orders=[];transforms={};end
old=find(orders==m,1);if ~isempty(old),T=transforms{old};return,end
N=zeros(m+1);D=zeros(1,m+1);
for k=0:m
    D(k+1)=nchoosek(m,k);
    for i=0:k,N(k+1,i+1)=nchoosek(k,i);end
end
require(all(N(:)<flintmax)&&all(D<flintmax),'Exact binomial integers exceed flintmax.');
T=sym(N)./sym(repmat(D,m+1,1));
orders(end+1)=m;transforms{end+1}=T;
end

function replayLeaves(B,paths,lower)
% Every internal node must have both halves in one and the same coordinate.
empty=cellfun(@isempty,paths);
if any(empty)
    require(numel(paths)==1,'Overlapping certificate leaves.');
    bound=rational(lower{1});
    require(isAlways(bound>0),'A retained leaf bound is not strictly positive.');
    verifyMinimum(B,bound);
    return
end
first=cellfun(@(p)p(1),paths);
if isequal(sort(unique(first(:))).','LR')
    axis=0; letters='LR';
elseif isequal(sort(unique(first(:))).','DU')
    axis=1; letters='DU';
else
    error('qstab:exact:IncompleteCoverage','Bernstein leaves do not cover the entire face.');
end
[lo,hi]=splitTensor(B,axis); children={lo,hi};
for k=1:2
    selected=first==letters(k);
    suffix=cellfun(@(p)p(2:end),paths(selected),'UniformOutput',false);
    replayLeaves(children{k},suffix,lower(selected));
end
end

function [lo,hi]=splitTensor(B,axis)
% Cached exact matrices implement dyadic de Casteljau subdivision. The
% matrices are universal basis transformations; no proof outcome is cached.
if axis==0
    [L,R]=splitTransform(size(B,1)-1);
    lo=L*B; hi=R*B;
else
    [L,R]=splitTransform(size(B,2)-1);
    lo=B*L.'; hi=B*R.';
end
end

function [L,R]=splitTransform(m)
persistent orders leftMatrices rightMatrices
if isempty(orders), orders=[]; leftMatrices={}; rightMatrices={}; end
old=find(orders==m,1);
if ~isempty(old), L=leftMatrices{old}; R=rightMatrices{old}; return; end
NL=zeros(m+1); NR=NL;
for k=0:m
    for i=0:k, NL(k+1,i+1)=nchoosek(k,i); end
    for i=k:m, NR(k+1,i+1)=nchoosek(m-k,i-k); end
end
DL=2.^(0:m).'; DR=2.^(m:-1:0).';
require(all(NL(:)<flintmax) && all(NR(:)<flintmax) && ...
    all(DL<flintmax) && all(DR<flintmax),'Exact split integers exceed flintmax.');
L=sym(NL)./sym(repmat(DL,1,m+1));
R=sym(NR)./sym(repmat(DR,1,m+1));
orders(end+1)=m; leftMatrices{end+1}=L; rightMatrices{end+1}=R;
end

function verifyMinimum(B,bound)
require(all(isAlways(B(:)>=bound)) && any(isAlways(B(:)==bound)), ...
    'Recomputed exact Bernstein minimum differs from the certificate.');
end

function q=rational(text)
text=char(text);
if isempty(regexp(text,'^-?[0-9]+(/[0-9]+)?$','once'))
    error('qstab:exact:Rational','Expected an integer or integer/integer rational literal.');
end
q=sym(text);
require(isempty(symvar(q)),'A coefficient is not an exact rational number.');
end

function tf=equalQ(a,b)
tf=isAlways(expand(a-b)==0);
end

function x=item(a,k)
if iscell(a), x=a{k}; else, x=a(k); end
end

function require(tf,message)
if ~all(tf(:)), error('qstab:exact:Verification','%s',message); end
end
