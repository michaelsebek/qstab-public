function R = zeroexcludeReport(data,varargin)
%ZEROEXCLUDEREPORT Internal implementation returning the full report.

ip=inputParser;
ip.FunctionName='qstab.zeroexclude';
addParameter(ip,'Uncertainty','',@localUncertaintyText);
addParameter(ip,'Reduction','none',@localText);
addParameter(ip,'MaxSubfamilies',1024,@(x)localInt(x,1));
addParameter(ip,'Variable','s',@localText);
addParameter(ip,'MaxCells',4096,@(x)localInt(x,1));
addParameter(ip,'MaxDepth',32,@(x)localInt(x,0));
addParameter(ip,'MaxTime',30,@(x)isnumeric(x)&&isscalar(x)&&isreal(x)&&x>=0&&~isnan(x));
addParameter(ip,'StoreCertificates',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Tolerance',1e-10,@(x)isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>0&&x<.01);
parse(ip,varargin{:}); o=ip.Results;
v=localVariable(o.Variable);
reduction=lower(strtrim(char(o.Reduction)));
if ~ismember(reduction,{'none','four-generator','kharitonov-cube','kharitonov-faces','quadratic-faces'})
    error('qstab:zeroexclude:Reduction', ...
        'Reduction must be none, four-generator, kharitonov-cube, kharitonov-faces, or quadratic-faces.');
end
F=qstab.familycomponents(data,'Uncertainty',o.Uncertainty);
if strcmp(reduction,'four-generator')
    R=qstab.internal.fourgenerator(F,o,v);
    return
elseif ismember(reduction,{'kharitonov-cube','kharitonov-faces','quadratic-faces'})
    R=qstab.internal.intervalreduce(F,o,v,reduction);
    return
end
% Joint parameters are retained for affine/one-quaternion families. For
% discrete coefficient boxes the same coefficient enters both E and O;
% their independent CT two-box shortcut must NOT be used on the unit circle.
if ismember(F.kind,{'affine','oneparameter'}) || ...
        (strcmp(F.kind,'interval') && ~strcmp(v,'s'))
    R=qstab.internal.affinezeroexclude(F,o,v);
    return
end
start=tic;
R=struct('status','undetermined','isStable',[],'method','normValueZeroExclusion', ...
  'uncertainty',F.kind,'variable',v,'testVariable',v,'degree',F.degree, ...
  'message','','assurance','floatingPointSufficientCertificate', ...
  'validated',false,'certificateComplete',false,'witness',[], ...
  'nominal',[],'frequencyLimit',[],'processedCells',0,'certifiedCells',0, ...
  'unresolvedCells',0,'maximumDepth',0,'elapsedTime',0, ...
  'unresolvedRegions',{{}},'certificates',{{}},'options',o);
if ~ismember(F.kind,{'interval','polytope'})
    R.message='Specialized zero exclusion currently supports interval and polytope families only.';
    R=localFinish(R,start,0,0,0,0);
    return
end
if strcmp(v,'zinv')
    R.testVariable='z';
end

% Change the time-variable convention if necessary, keeping the declared grade.
if strcmp(F.kind,'polytope')
    V=F.vertices;
    if strcmp(v,'zinv'), V=V(:,end:-1:1,:); end
    for j=1:size(V,3)
        if all(V(:,end,j)==0)
            R.status='assumptionFailed'; R.assurance='floatingPointAssumptionWitness';
            R.witness=struct('components',V(:,:,j),'originalWeights',localBasis(j,size(V,3)), ...
                'kind','zeroLeadingCoefficient');
            R.message='An admissible generator has zero leading coefficient in the tested variable.';
            R=localFinish(R,start,0,0,0,0);
            return
        end
    end
    % Certify the FULL leading hull before using zero-exclusion continuity.
    % In particular, nonzero endpoint coefficients alone are insufficient.
    R.leading=qstab.internal.leadbound(F,v,o.Tolerance,256, ...
        max(0,o.MaxTime-toc(start)));
    if strcmp(R.leading.status,'zeroWitness')
        R.status='assumptionFailed';R.assurance='floatingPointAssumptionWitness';
        R.witness=R.leading.witness;R.message=R.leading.message;
        R=localFinish(R,start,0,0,0,0);return
    elseif ~R.leading.isSeparated
        R.message=['Full leading hull is unresolved: ',R.leading.message];
        R.unresolvedRegions={struct('kind','leadingCoefficientUnresolved', ...
            'vertexWeights',eye(size(V,3)),'reason',R.message)};
        R=localFinish(R,start,0,0,1,0);return
    end
    mlead=R.leading.lowerBound;
    allUpper=zeros(1,F.degree);
    for j=1:F.degree
        for k=1:size(V,3)
            allUpper(j)=max(allUpper(j),norm(V(:,j,k)));
        end
    end
    nodes={struct('lo',0,'hi',0,'vertices',V,'weights',eye(size(V,3)),'depth',0)};
else
    L=F.intervalLower;
    U=F.intervalUpper;
    leadN=F.degree+1;
    minAbs=max(0,max(L(:,leadN),-U(:,leadN)));
    mlead=norm(minAbs);
    if mlead==0
        R.status='assumptionFailed'; R.assurance='floatingPointAssumptionWitness';
        W=(L+U)/2;
        W(:,leadN)=min(max(zeros(4,1),L(:,leadN)),U(:,leadN));
        R.witness=struct('components',W,'kind','zeroLeadingCoefficient');
        R.message='The Cartesian interval contains a zero leading quaternion.';
        R=localFinish(R,start,0,0,0,0);
        return
    end
    allUpper=zeros(1,F.degree);
    for j=1:F.degree
        allUpper(j)=norm(max(abs(L(:,j)),abs(U(:,j))));
    end
    nodes={struct('lo',0,'hi',0,'elo',-inf(4,1),'ehi',inf(4,1), ...
        'olo',-inf(4,1),'ohi',inf(4,1),'depth',0)};
end

if F.degree==0
    R.status='stable'; R.isStable=true; R.certificateComplete=true;
    R.message='Every member is a nonzero constant polynomial.';
    R=localFinish(R,start,0,1,0,0);
    return
end

if ~(mlead>0) || ~isfinite(mlead)
    R.message='Uniform leading-coefficient bound is unavailable.';
    R=localFinish(R,start,0,0,0,0);
    return
end
if strcmp(v,'s')
    freqMax=1+max(allUpper)/mlead;
else
    freqMax=pi;  % upper unit semicircle, z=exp(i*angle)
end
if ~isfinite(freqMax)
    R.message='Uniform frequency cutoff exceeds representable range.';
    R=localFinish(R,start,0,0,0,0);
    return
end
R.frequencyLimit=freqMax;
nodes{1}.hi=freqMax;

% Check a real admissible member before invoking zero-exclusion continuity.
if strcmp(F.kind,'polytope')
    Cp=sum(V/size(V,3),3);
    weights=ones(size(V,3),1)/size(V,3);
else
    Cp=(L+U)/2;
    weights=[];
end
R.nominal=localNominal(Cp,R.testVariable,o.Tolerance);
if strcmp(R.nominal.status,'notStable')
    R.status='notStable'; R.isStable=false;
    R.assurance='floatingPointCounterexample';
    R.witness=struct('components',Cp,'originalWeights',weights, ...
        'kind','actualAdmissiblePolynomial','nominal',R.nominal);
    R.message='An admissible central polynomial is not strictly stable.';
    R=localFinish(R,start,0,0,0,0);
    return
elseif ~strcmp(R.nominal.status,'stable')
    R.message='No strictly stable reference was verified by the nominal test.';
    R=localFinish(R,start,0,0,0,0);
    return
end

processed=0; certified=0; depthMax=0; lastReason='';
while ~isempty(nodes)
    if processed>=o.MaxCells
        lastReason='MaxCells reached'; break
    end
    if toc(start)>=o.MaxTime
        lastReason='MaxTime reached'; break
    end
    N=nodes{end}; nodes(end)=[];
    processed=processed+1;
    depthMax=max(depthMax,N.depth);
    if strcmp(F.kind,'polytope')
        [excluded,diagnostic]=localPolyCell(N,R.testVariable,o.Tolerance);
    else
        [excluded,diagnostic,N]=localIntervalCell(N,L,U,o.Tolerance);
    end
    if excluded
        certified=certified+1;
        if o.StoreCertificates
            R.certificates{end+1,1}=struct('cell',localRegion(N,F.kind, ...
                'Accepted zero-exclusion region'),'diagnostic',diagnostic); %#ok<AGROW>
        end
        continue
    end
    lastReason=diagnostic.reason;
    if strcmp(F.kind,'polytope')
        cp=sum(N.vertices/size(N.vertices,3),3);
        cweights=mean(N.weights,2);
        candidate=localNominal(cp,R.testVariable,o.Tolerance);
        if strcmp(candidate.status,'notStable')
            R.status='notStable'; R.isStable=false;
            R.assurance='floatingPointCounterexample';
            R.witness=struct('components',cp,'originalWeights',cweights, ...
                'kind','actualAdmissiblePolynomial','nominal',candidate);
            R.message='A polynomial from an admissible parameter simplex is not stable.';
            R=localPending(R,nodes,'Not checked after a counterexample');
            R=localFinish(R,start,processed,certified,numel(R.unresolvedRegions),depthMax);
            return
        end
    end
    if N.depth>=o.MaxDepth
        R.unresolvedRegions{end+1,1}=localRegion(N,F.kind,'MaxDepth reached'); %#ok<AGROW>
        continue
    end
    if strcmp(F.kind,'polytope')
        if diagnostic.centerDistance>0 && ~diagnostic.needsParameterSplit && ...
                diagnostic.freqVariation >= 0.5*diagnostic.centerDistance
            [A,B,ok]=localFreqSplit(N);
        else
            [A,B,ok]=localPolySplit(N);
            if ~ok, [A,B,ok]=localFreqSplit(N); end
        end
    else
        if diagnostic.freqExcess>0.4*max(diagnostic.totalSpan,realmin)
            [A,B,ok]=localFreqSplit(N);
        else
            [A,B,ok]=localIntervalSplit(N);
            if ~ok, [A,B,ok]=localFreqSplit(N); end
        end
    end
    if ~ok
        R.unresolvedRegions{end+1,1}=localRegion(N,F.kind,'Cannot further subdivide'); %#ok<AGROW>
        continue
    end
    nodes{end+1}=B; %#ok<AGROW>
    nodes{end+1}=A; %#ok<AGROW>
end
R=localPending(R,nodes,lastReason);
if isempty(R.unresolvedRegions) && isempty(nodes)
    R.status='stable'; R.isStable=true; R.certificateComplete=true;
    R.message='Norm-polynomial zero was excluded on every covering frequency/parameter cell.';
else
    R.status='undetermined'; R.isStable=[];
    R.assurance='floatingPointUndetermined';
    R.message=['No complete zero-exclusion certificate within the budget: ',lastReason];
end
R=localFinish(R,start,processed,certified,numel(R.unresolvedRegions),depthMax);
end

function [excluded,D]=localPolyCell(N,variable,tol)
% The degree-two Bernstein controls for N_{sum alpha_j p_j}.
V=N.vertices; m=size(V,3); n=size(V,2)-1;
controls=complex(zeros(1,m*(m+1)/2));
controlsLo=controls; controlsHi=controls;
Lmax=0;
center=(N.lo+N.hi)/2;
if strcmp(variable,'s')
    z=1i*center;
    R=max(1,N.hi);
else
    z=exp(1i*center);
    R=1;
end
kindex=0;
for a=1:m
    for b=a:m
        kindex=kindex+1;
        beta=zeros(1,2*n+1);
        for nu=1:4
            beta=beta+conv(V(nu,:,a),V(nu,:,b));
        end
        if any(~isfinite(beta))
            excluded=false;
            D=struct('reason','Polarization coefficients overflowed', ...
               'centerDistance',0,'freqVariation',inf,'needsParameterSplit',false);
            return
        end
        q=0;
        for ell=numel(beta):-1:1
            q=q*z+beta(ell);
        end
        controls(kindex)=q;
        if strcmp(variable,'s')
            zlo=1i*N.lo; zhi=1i*N.hi;
        else
            zlo=exp(1i*N.lo); zhi=exp(1i*N.hi);
        end
        controlsLo(kindex)=localAscEval(beta,zlo);
        controlsHi(kindex)=localAscEval(beta,zhi);
        powers=R.^(0:(2*n-1));
        if any(~isfinite(powers))
            Lmax=inf;
        else
            Lmax=max(Lmax,sum((1:2*n).*abs(beta(2:end)).*powers));
        end
    end
end
[dist,wellSeparated]=localHullDistance(controls,tol);
[distLo,~]=localHullDistance(controlsLo,tol);
[distHi,~]=localHullDistance(controlsHi,tol);
freqError=Lmax*(N.hi-N.lo)/2;
scale=max([1,abs(controls)]);
errorAllowance=tol*(1+scale+freqError);
excluded=wellSeparated && isfinite(freqError) && ...
    dist > freqError+errorAllowance;
D=struct('reason','Bernstein hull/frequency enclosure still meets zero', ...
  'centerDistance',dist,'freqVariation',freqError,'valueScale',scale, ...
  'needsParameterSplit',distLo==0 || distHi==0);
if excluded
    D.reason='Strict Bernstein hull/frequency separation';
end
end

function q=localAscEval(beta,z)
q=0;
for j=numel(beta):-1:1
    q=q*z+beta(j);
end
end

function [dist,wellSeparated]=localHullDistance(v,tol)
% Exact Euclidean distance to conv(v) in 2D, except for floating roundoff.
v=v(:);
dist=0; wellSeparated=false;
if any(~isfinite(real(v))|~isfinite(imag(v))) || isempty(v), return; end
if any(abs(v)<tol*max(1,max(abs(v)))), return; end
angles=sort(mod(angle(v),2*pi));
gaps=diff([angles; angles(1)+2*pi]);
% A strict open semicircle exists iff the maximum circular gap exceeds pi.
if max(gaps)<=pi+32*eps*pi, return; end
wellSeparated=true;
dist=min(abs(v));
for i=1:numel(v)
    for j=i+1:numel(v)
        d=v(j)-v(i);
        length2=abs(d)^2;
        if length2==0, continue; end
        alpha=max(0,min(1,-real(conj(d)*v(i))/length2));
        dist=min(dist,abs(v(i)+alpha*d));
    end
end
end

function [excluded,D,N]=localIntervalCell(N,L,U,tol)
% Continuous boundary w >= 0: independent exact E/O boxes at each w.
[el,eu,ol,ou]=localEvenOddBounds(L,U,N.lo,N.hi);
[eml,emu,oml,omu]=localEvenOddBounds(L,U,(N.lo+N.hi)/2,(N.lo+N.hi)/2);
N.elo=max(N.elo,el); N.ehi=min(N.ehi,eu);
N.olo=max(N.olo,ol); N.ohi=min(N.ohi,ou);
if any(N.elo>N.ehi) || any(N.olo>N.ohi)
    excluded=true;
    D=struct('reason','Empty evaluation-box intersection', ...
       'freqExcess',0,'totalSpan',0);
    return
end
[eMin,eMax]=localSquaredRange(N.elo,N.ehi);
[oMin,oMax]=localSquaredRange(N.olo,N.ohi);
realLo=eMin-oMax; realHi=eMax-oMin;
prodLo=0; prodHi=0;
for j=1:4
    prod=[N.elo(j)*N.olo(j),N.elo(j)*N.ohi(j), ...
          N.ehi(j)*N.olo(j),N.ehi(j)*N.ohi(j)];
    prodLo=prodLo+min(prod); prodHi=prodHi+max(prod);
end
imagLo=2*prodLo; imagHi=2*prodHi;
margin=tol*(1+max([abs(realLo),abs(realHi),abs(imagLo),abs(imagHi)]));
excluded=isfinite(margin) && ...
  (realLo>margin || realHi< -margin || imagLo>margin || imagHi< -margin);
span=sum(N.ehi-N.elo)+sum(N.ohi-N.olo);
midSpan=sum(max(0,min(N.ehi,emu)-max(N.elo,eml)))+ ...
        sum(max(0,min(N.ohi,omu)-max(N.olo,oml)));
D=struct('reason','Even/odd interval norm image still meets zero', ...
 'freqExcess',max(0,span-midSpan),'totalSpan',span, ...
 'realInterval',[realLo realHi],'imagInterval',[imagLo imagHi]);
if excluded
    D.reason='Strict real or imaginary norm-value interval separation';
end
end

function [el,eu,ol,ou]=localEvenOddBounds(L,U,w0,w1)
el=zeros(4,1); eu=zeros(4,1); ol=zeros(4,1); ou=zeros(4,1);
for k=0:size(L,2)-1
    a=w0^k; b=w1^k;
    if mod(floor(k/2),2)==1
        a=-a; b=-b;
    end
    if k==0, a=1; b=1; end
    for nu=1:4
        v=[a*L(nu,k+1),a*U(nu,k+1),b*L(nu,k+1),b*U(nu,k+1)];
        if mod(k,2)==0
            el(nu)=el(nu)+min(v);
            eu(nu)=eu(nu)+max(v);
        else
            ol(nu)=ol(nu)+min(v);
            ou(nu)=ou(nu)+max(v);
        end
    end
end
end

function [minVal,maxVal]=localSquaredRange(l,u)
minAbs=max(0,max(l,-u));
minVal=sum(minAbs.^2);
maxVal=sum(max(abs(l),abs(u)).^2);
end

function [A,B,ok]=localFreqSplit(N)
mid=N.lo+(N.hi-N.lo)/2;
ok=mid>N.lo && mid<N.hi && isfinite(mid);
A=N; B=N;
if ok
    A.hi=mid; B.lo=mid;
    A.depth=N.depth+1; B.depth=N.depth+1;
end
end

function [A,B,ok]=localPolySplit(N)
V=N.vertices; m=size(V,3);
A=N; B=N;
longest=0; iBest=0; jBest=0;
for i=1:m
    for j=i+1:m
        d=norm(V(:,:,i)-V(:,:,j),'fro');
        if d>longest, longest=d; iBest=i; jBest=j; end
    end
end
ok=longest>0 && isfinite(longest);
if ~ok, return; end
midV=0.5*V(:,:,iBest)+0.5*V(:,:,jBest);
midW=0.5*N.weights(:,iBest)+0.5*N.weights(:,jBest);
if isequal(midV,V(:,:,iBest)) || isequal(midV,V(:,:,jBest))
    ok=false; return
end
A.vertices(:,:,iBest)=midV;
B.vertices(:,:,jBest)=midV;
A.weights(:,iBest)=midW;
B.weights(:,jBest)=midW;
A.depth=N.depth+1; B.depth=N.depth+1;
end

function [A,B,ok]=localIntervalSplit(N)
A=N; B=N;
width=[N.ehi-N.elo; N.ohi-N.olo];
[wid,index]=max(width);
ok=isfinite(wid) && wid>0;
if ~ok, return; end
if index<=4
    l=N.elo(index); u=N.ehi(index);
else
    l=N.olo(index-4); u=N.ohi(index-4);
end
mid=l+(u-l)/2;
ok=mid>l && mid<u;
if ~ok, return; end
if index<=4
    A.ehi(index)=mid; B.elo(index)=mid;
else
    A.ohi(index-4)=mid; B.olo(index-4)=mid;
end
A.depth=N.depth+1; B.depth=N.depth+1;
end

function out=localRegion(N,kind,reason)
out=struct('kind',kind,'frequencyBounds',[N.lo,N.hi], ...
 'depth',N.depth,'reason',reason,'vertexWeights',[], ...
 'evenBox',[],'oddBox',[]);
if strcmp(kind,'polytope')
    out.vertexWeights=N.weights;
else
    out.evenBox=[N.elo,N.ehi];
    out.oddBox=[N.olo,N.ohi];
end
end

function R=localPending(R,nodes,reason)
for k=1:numel(nodes)
    R.unresolvedRegions{end+1,1}=localRegion(nodes{k},R.uncertainty,reason); %#ok<AGROW>
end
end

function R=localFinish(R,started,processed,certified,unresolved,maxDepth)
if strcmp(R.status,'undetermined')
    R.assurance='floatingPointUndetermined';
end
R.processedCells=processed;
R.certifiedCells=certified;
R.unresolvedCells=unresolved;
R.maximumDepth=maxDepth;
R.elapsedTime=toc(started);
end

function N=localNominal(C,variable,tol)
S=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
try
    N=qstab.internal.isstableReport(S,'Variable',variable,'Tolerance',tol);
catch ME
    N=struct('status','undetermined','message',ME.message);
end
end

function w=localBasis(j,m)
w=zeros(m,1); w(j)=1;
end

function tf=localUncertaintyText(x)
% Empty character text is the documented sentinel for family inference.
% robuststable explicitly forwards its default ''. It must be accepted in
% exactly the same way as an omitted Uncertainty option; [] is not text.
tf=localText(x)||(ischar(x)&&isempty(x));
end

function tf=localText(x)
tf=(ischar(x)&&isrow(x))||(isstring(x)&&isscalar(x));
end
function tf=localInt(x,lowerLimit)
tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x==floor(x)&&x>=lowerLimit;
end
function v=localVariable(x)
v=lower(strtrim(char(x)));
if ismember(v,{'continuous','hurwitz'}),v='s';end
if ismember(v,{'discrete','schur'}),v='z';end
if ismember(v,{'delay','z^-1','z^{-1}','d'}),v='zinv';end
if ~ismember(v,{'s','z','zinv'})
    error('qstab:zeroexclude:Variable','Variable must be s, z, or zinv.');
end
end
