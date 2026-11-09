function R = affinezeroexclude(F,o,variable)
%AFFINEZEROEXCLUDE Joint-parameter norm-image enclosures for affine boxes.
% All coefficients retain the SAME real parameter vector. At a frequency,
% write p=pc+sum xi_j*hj, -1<=xi_j<=1, and bound the exact quadratic
% G=N_pc+2*sum xi_j B(pc,hj)+sum xi_j^2 N_hj+
%       2*sum_(j<k) xi_j xi_k B(hj,hk).
% An outer rectangle plus a whole-frequency derivative disk is used only
% for sufficient exclusion. Its intersection with zero is NOT a witness.
started=tic;tv=variable;if strcmp(variable,'zinv'),tv='z';end
R=struct('status','undetermined','isStable',[], ...
    'method','affineNormValueZeroExclusion','uncertainty',F.kind, ...
    'variable',variable,'testVariable',tv,'degree',F.degree, ...
    'message','','assurance','floatingPointSufficientCertificate', ...
    'validated',false,'certificateComplete',false,'witness',[], ...
    'nominal',[],'frequencyLimit',[],'processedCells',0,'certifiedCells',0, ...
    'unresolvedCells',0,'maximumDepth',0,'elapsedTime',0, ...
    'unresolvedRegions',{{}},'certificates',{{}},'options',o);
R.parameterNames=F.parameterNames;R.parameterIndex=F.parameterIndex;
R.leading=qstab.internal.leadbound(F,variable,o.Tolerance,256, ...
    max(0,o.MaxTime-toc(started)));
if strcmp(R.leading.status,'zeroWitness')
    R.status='assumptionFailed';R.assurance='floatingPointAssumptionWitness';
    R.witness=R.leading.witness;R.message=R.leading.message;
    R=localFinish(R,started);return
elseif ~R.leading.isSeparated
    R.message=['Invariant degree was not established: ',R.leading.message];
    R.unresolvedRegions={struct('kind','leadingCoefficientUnresolved', ...
        'parameterBounds',F.bounds,'frequencyBounds',[], ...
        'reason',R.message)};
    R=localFinish(R,started);return
end
if F.degree==0
    R.status='stable';R.isStable=true;R.certificateComplete=true;
    R.message='Every admissible polynomial is a nonzero constant.';
    R=localFinish(R,started);return
end
B=F.base;G=F.generators;
if strcmp(variable,'zinv'),B=B(:,end:-1:1);G=G(:,end:-1:1,:);end
mid=.5*F.bounds(:,1)+.5*F.bounds(:,2);
rad=.5*F.bounds(:,2)-.5*F.bounds(:,1);
C=localValue(B,G,mid);upper=abs(C);
for j=1:numel(mid),upper=upper+rad(j)*abs(G(:,:,j));end
if any(~isfinite(upper(:)))
    R.message='The uniform coefficient enclosure overflowed.';
    R=localWholeUnresolved(R,F,R.message);R=localFinish(R,started);return
end
if strcmp(tv,'s')
    bound=zeros(1,F.degree);
    for k=1:F.degree,bound(k)=norm(upper(:,k));end
    limit=1+max(bound)/R.leading.lowerBound;
else
    limit=pi;
end
if ~isfinite(limit)
    R.message='The full-family root bound is outside the representable range.';
    R=localWholeUnresolved(R,F,R.message);R=localFinish(R,started);return
end
R.frequencyLimit=limit;
R.nominal=localNominal(C,tv,o.Tolerance);
if strcmp(R.nominal.status,'notStable')
    R.status='notStable';R.isStable=false;R.assurance='floatingPointCounterexample';
    R.witness=localWitness(F,mid,variable,R.nominal);
    R.message='An admissible parameter midpoint is not stable.';
    R=localFinish(R,started);return
elseif ~strcmp(R.nominal.status,'stable')
    R.message='The reference midpoint is not resolved as strictly stable.';
    R=localWholeUnresolved(R,F,R.message);R=localFinish(R,started);return
end
stack={struct('bounds',F.bounds,'lo',0,'hi',limit,'depth',0,'nominal',R.nominal)};
lastReason='';
while ~isempty(stack)
    if R.processedCells>=o.MaxCells,lastReason='MaxCells reached';break;end
    if toc(started)>=o.MaxTime,lastReason='MaxTime reached';break;end
    N=stack{end};stack(end)=[];
    R.processedCells=R.processedCells+1;R.maximumDepth=max(R.maximumDepth,N.depth);
    [excluded,D]=localEnclosure(B,G,N,tv,o.Tolerance);
    if excluded
        R.certifiedCells=R.certifiedCells+1;
        if o.StoreCertificates
            R.certificates{end+1,1}=struct('cell',localRegion(N,'Accepted affine cell'), ...
                'diagnostic',D); %#ok<AGROW>
        end
        continue
    end
    theta=.5*N.bounds(:,1)+.5*N.bounds(:,2);
    if isempty(N.nominal),N.nominal=localNominal(localValue(B,G,theta),tv,o.Tolerance);end
    if strcmp(N.nominal.status,'notStable')
        R.status='notStable';R.isStable=false;R.assurance='floatingPointCounterexample';
        R.witness=localWitness(F,theta,variable,N.nominal);
        R.message='An admissible correlated parameter member is not stable.';
        R=localPending(R,stack,'Not evaluated after a counterexample');
        R=localFinish(R,started);return
    end
    if N.depth>=o.MaxDepth
        R.unresolvedRegions{end+1,1}=localRegion(N,'MaxDepth reached'); %#ok<AGROW>
        lastReason='MaxDepth reached';continue
    end
    [A,C,ok]=localSplit(N,G,D);
    if ~ok
        R.unresolvedRegions{end+1,1}=localRegion(N,D.reason); %#ok<AGROW>
        lastReason=D.reason;continue
    end
    stack{end+1}=C;stack{end+1}=A; %#ok<AGROW>
end
R=localPending(R,stack,lastReason);
if isempty(R.unresolvedRegions)
    R.status='stable';R.isStable=true;R.certificateComplete=true;
    R.message='Every joint-parameter/frequency cell strictly excludes a norm-polynomial zero.';
else
    R.assurance='floatingPointUndetermined';
    R.message=['Joint-parameter zero exclusion is incomplete: ',lastReason];
end
R=localFinish(R,started);
end

function [excluded,D]=localEnclosure(B,G,N,variable,tol)
t=.5*N.bounds(:,1)+.5*N.bounds(:,2);
r=.5*N.bounds(:,2)-.5*N.bounds(:,1);
C=localValue(B,G,t);m=numel(t);n=size(C,2)-1;
f=.5*N.lo+.5*N.hi;
if strcmp(variable,'s'),z=1i*f;rho=max(1,N.hi);else,z=exp(1i*f);rho=1;end
H=bsxfun(@times,G,reshape(r,1,1,[]));
y=localEval(C,z);Dval=complex(zeros(4,m));
for j=1:m,Dval(:,j)=localEval(H(:,:,j),z);end
q=y.'*y;linear=2*(Dval.'*y);Q=Dval.'*Dval;diagonal=diag(Q);
off=triu(Q,1);
lo=[real(q)-sum(abs(real(linear)))+sum(min(0,real(diagonal)))-2*sum(abs(real(off(:)))); ...
    imag(q)-sum(abs(imag(linear)))+sum(min(0,imag(diagonal)))-2*sum(abs(imag(off(:))))];
hi=[real(q)+sum(abs(real(linear)))+sum(max(0,real(diagonal)))+2*sum(abs(real(off(:)))); ...
    imag(q)+sum(abs(imag(linear)))+sum(max(0,imag(diagonal)))+2*sum(abs(imag(off(:))))];
U=abs(C)+sum(abs(H),3);beta=zeros(1,2*n+1);
for nu=1:4,beta=beta+conv(U(nu,:),U(nu,:));end
L=sum((1:2*n).*beta(2:end).*rho.^(0:2*n-1));
variation=L*(N.hi-N.lo)/2;
near=max(0,max(lo,-hi));distance=norm(near);
scale=max([1;abs(lo);abs(hi);variation]);
allowance=(tol+256*eps)*scale;
finite=all(isfinite([lo;hi;variation;distance;allowance]));
excluded=finite && distance>variation+allowance;
D=struct('realRange',[lo(1),hi(1)],'imagRange',[lo(2),hi(2)], ...
    'centerDistance',distance,'freqVariation',variation, ...
    'parameterSpan',max(hi-lo)/2,'errorAllowance',allowance, ...
    'gap',distance-variation-allowance,'reason','Joint quadratic enclosure still meets zero');
if ~finite,D.reason='Nonfinite joint-parameter enclosure';end
if excluded,D.reason='Strict joint-parameter rectangle/frequency separation';end
end

function y=localEval(C,z)
y=complex(zeros(size(C,1),1));
for k=size(C,2):-1:1,y=y*z+C(:,k);end
end
function C=localValue(B,G,t)
C=B;for j=1:numel(t),C=C+t(j)*G(:,:,j);end
end
function S=localPolynomial(C)
S=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
function S=localNominal(C,variable,tol)
if any(~isfinite(C(:)))
    S=struct('status','undetermined','isStable',[]);return
end
S=qstab.internal.isstableReport(localPolynomial(C),'Variable',variable,'Tolerance',tol);
end
function W=localWitness(F,t,variable,nominal)
C=localValue(F.base,F.generators,t);Ct=C;tv=variable;
if strcmp(variable,'zinv'),Ct=C(:,end:-1:1);tv='z';end
W=struct('kind','actualAdmissiblePolynomial','components',Ct, ...
    'originalComponents',C,'parameter',t,'originalWeights',[], ...
    'nominal',nominal,'coefficientVariable',tv, ...
    'originalCoefficientVariable',variable);
end
function r=localRegion(N,reason)
r=struct('kind','affine','parameterBounds',N.bounds, ...
    'frequencyBounds',[N.lo,N.hi],'depth',N.depth,'reason',reason);
end
function R=localWholeUnresolved(R,F,reason)
R.unresolvedRegions={struct('kind','affine','parameterBounds',F.bounds, ...
    'frequencyBounds',R.frequencyLimit,'depth',0,'reason',reason)};
end
function R=localPending(R,nodes,reason)
for j=1:numel(nodes),R.unresolvedRegions{end+1,1}=localRegion(nodes{j},reason);end
end
function R=localFinish(R,started)
R.unresolvedCells=numel(R.unresolvedRegions);R.elapsedTime=toc(started);
end
function [A,B,ok]=localSplit(N,G,D)
A=N;B=N;ok=false;
mid=.5*N.lo+.5*N.hi;canFreq=mid>N.lo && mid<N.hi;
width=N.bounds(:,2)-N.bounds(:,1);scores=zeros(numel(width),1);
for j=1:numel(width),scores(j)=width(j)*norm(reshape(G(:,:,j),[],1));end
[score,j]=max(scores);canParam=~isempty(score)&&score>0;
if canFreq && (~canParam || D.freqVariation>=D.parameterSpan/2)
    A.hi=mid;B.lo=mid;ok=true;
elseif canParam
    t=.5*N.bounds(j,1)+.5*N.bounds(j,2);
    if t>N.bounds(j,1) && t<N.bounds(j,2)
        A.bounds(j,2)=t;B.bounds(j,1)=t;A.nominal=[];B.nominal=[];ok=true;
    end
elseif canFreq
    A.hi=mid;B.lo=mid;ok=true;
end
if ok,A.depth=N.depth+1;B.depth=N.depth+1;end
end
