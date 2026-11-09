function R = discreteray(X0,X1,u,o)
%DISCRETERAY Connected open Schur interval along t=alpha*u.
% An unnormalized companion descriptor preserves affine leading coefficients.
% Its real quaternion-Hermitian generalized Stein operator is quadratic in
% the real ray parameter. QZ supplies candidates; ambiguous earlier contacts
% block later endpoints. Degree/regularity loss is an additional event.
% No frequency or ray grid is used as a certificate; no epsilon pivot fix.
typ=lower(char(o.Type));side=lower(char(o.Side));variable=lower(char(o.Variable));
if ~any(strcmp(typ,{'polynomial','matrix'}))||~any(strcmp(side,{'left','right'}))
    error('qstab:raymargin:Option','Invalid Type or Side.');
end
if ~(isnumeric(o.Tolerance)&&isreal(o.Tolerance)&&isscalar(o.Tolerance)&&isfinite(o.Tolerance)&&o.Tolerance>0&&o.Tolerance<1)
    error('qstab:raymargin:Tolerance','Tolerance must lie strictly between zero and one.');
end
if ~(isscalar(o.VerifyBoundary)&&(islogical(o.VerifyBoundary)|| ...
        (isnumeric(o.VerifyBoundary)&&isreal(o.VerifyBoundary)&&any(o.VerifyBoundary==[0 1]))))
    error('qstab:raymargin:VerifyBoundary','VerifyBoundary must be a logical scalar.');
end
if ~integer(o.MaxGuardCells,0)||~integer(o.MaxGuardDepth,0)|| ...
        ~(isnumeric(o.MaxGuardTime)&&isreal(o.MaxGuardTime)&&isscalar(o.MaxGuardTime)&&~isnan(o.MaxGuardTime)&&o.MaxGuardTime>=0)
    error('qstab:raymargin:GuardBudget','Invalid whole-segment guard budget.');
end
[a,b,c,d]=qstab.components(u);if numel(a)~=1,error('qstab:raymargin:Direction','u must be scalar.');end
uv=[a;b;c;d];un=norm(uv);
if ~(un>0&&isfinite(un)),error('qstab:raymargin:Direction','u must be finite and nonzero.');end
unit=uv/un;
R=struct('schema','qstab.discrete-ray.1','status','undetermined', ...
    'method','quaternion-hermitian-generalized-Stein-quadratic', ...
    'assurance','floating-point','validated',false,'type',typ,'side',side, ...
    'variable',variable,'testVariable','z','directionComponents',uv,'directionNorm',un, ...
    'interval',[NaN NaN],'candidateInterval',[NaN NaN], ...
    'negative',NaN,'positive',NaN,'negativeStatus','undetermined','positiveStatus','undetermined', ...
    'positiveParameterNorm',NaN,'negativeParameterNorm',NaN, ...
    'negativeBracket',[NaN NaN],'positiveBracket',[NaN NaN], ...
    'negativeEvent','none','positiveEvent','none','negativeWitness',[],'positiveWitness',[], ...
    'message','','diagnostics',struct());
M=[];leadCandidate=NaN;leadExact=false;lead0=[];lead1=[];
if strcmp(typ,'polynomial')
    M=qstab.internal.discretepair(X0,X1,variable,o.Grade);n=M.grade;R.grade=n;
    nominal=qstab.internal.isstableReport(polyParts(M.C0),'Variable','z','Grade',n,'Tolerance',o.Tolerance);
    R.nominal=nominal;
    if ~strcmp(nominal.status,'stable')
        if any(strcmp(nominal.status,{'notStable','assumptionFailed'})),R.status='assumptionFailed';end
        R.message='Require a verified Schur nominal polynomial at the declared forward grade.';return
    end
    if ~M.hasPerturbation,R=unbounded(R,'The polynomial perturbation is identically zero.');return;end
    B=zeros(size(M.S1));
    for j=1:n+1
        if strcmp(side,'left'),B(:,j)=mul(unit,M.S1(:,j));else,B(:,j)=mul(M.S1(:,j),unit);end
    end
    logParameterScale=M.logScaleRatio-log(un);C=M.S0;
    % Exact structural proportionality is not inferred from a tolerance.
    [~,pivot]=max(abs(C(:)));factor=B(pivot)/C(pivot);
    exactSame=isequal(M.C1,M.C0)||isequal(M.C1,-M.C0);
    if all(uv(2:4)==0)&&exactSame&&factor~=0&&all(B(:)==factor*C(:))
        betaLoss=-1/factor;alphaLoss=scaleval(betaLoss,logParameterScale);
        if ~isfinite(alphaLoss)||alphaLoss==0
            R.message='The scalar regularity-loss parameter is outside floating-point range.';return
        end
        R=unbounded(R,'The family is a real scalar multiple of the nominal polynomial until its zero-polynomial event.');
        [W,ok]=makeWitness(X0,X1,M,uv,alphaLoss,side,typ,'degreeLoss',o.Tolerance);
        if ~ok,R.status='undetermined';R.interval=[NaN NaN];R.positive=NaN;R.negative=NaN; ...
                R.positiveStatus='undetermined';R.negativeStatus='undetermined';return;end
        if alphaLoss>0,R.positive=alphaLoss;R.interval(2)=alphaLoss;R.positiveEvent='degreeLoss';R.positiveWitness=W;R.positiveParameterNorm=norm(W.parameterComponents);
        else,R.negative=alphaLoss;R.interval(1)=alphaLoss;R.negativeEvent='degreeLoss';R.negativeWitness=W;R.negativeParameterNorm=norm(W.parameterComponents);end
        R.candidateInterval=R.interval;R.diagnostics.structuralProportionality=true;return
    end
    if n==0
        % A nonzero quaternion scalar has no finite modes. Only exact zero of
        % its affine real ray matters; near-collinearity is not rounded away.
        [leadCandidate,leadExact,leadNear]=leadingContact(C(:,end),B(:,end));
        if ~isfinite(leadCandidate),R=unbounded(R,'A nonzero affine quaternion constant has no real zero on this ray.');return;end
        if ~leadExact&&leadNear,R.message='A near-collinear constant-family zero is numerically ambiguous.';return;end
        R.message='The constant-family regularity event is unresolved.';return
    end
    if n>16,error('qstab:raymargin:OrderLimit','Discrete dense Stein rays currently support grade/order at most 16.');end
    [A0,A1,E0,E1]=companionPair(C,B,n);
    [leadCandidate,leadExact,~]=leadingContact(C(:,end),B(:,end));
    lead0=C(:,end);lead1=B(:,end);
else
    if ~strcmpi(o.Domain,'discrete')||~any(strcmp(variable,{'s','z'}))||~isempty(o.Grade)
        error('qstab:raymargin:Domain','Matrices use Domain=discrete; zinv and polynomial Grade do not describe a state-update matrix.');
    end
    R.variable='z';
    [a,b,c,d]=qstab.components(X0);Q0=cat(3,a,b,c,d);
    [a,b,c,d]=qstab.components(X1);Q1=cat(3,a,b,c,d);
    if ndims(a)>2||isempty(a)||size(a,1)~=size(a,2)||~isequal(size(Q0),size(Q1))
        error('qstab:raymargin:MatrixShape','A0 and A1 must be square matrices of equal size.');
    end
    n=size(a,1);if n>16,error('qstab:raymargin:OrderLimit','Discrete dense Stein rays currently support grade/order at most 16.');end
    nominal=qstab.internal.isstableReport(X0,'Type','matrix','Domain','discrete','Tolerance',o.Tolerance);R.nominal=nominal;
    if ~strcmp(nominal.status,'stable')
        if any(strcmp(nominal.status,{'notStable','assumptionFailed'})),R.status='assumptionFailed';end
        R.message='Require a verified Schur nominal state matrix.';return
    end
    slopeScale=max(abs(Q1(:)));
    if slopeScale==0,R=unbounded(R,'The matrix perturbation is identically zero.');return;end
    B=zeros(size(Q1));
    for i=1:n,for j=1:n
        q=reshape(Q1(i,j,:),4,1)/slopeScale;
        if strcmp(side,'left'),q=mul(unit,q);else,q=mul(q,unit);end
        B(i,j,:)=reshape(q,1,1,4);
    end;end
    if fixedTriangular(Q0,B)
        R=unbounded(R,'A common triangular family has an unchanged quaternion diagonal and hence unchanged right eigenvalue classes.');return
    end
    logParameterScale=-log(slopeScale)-log(un);
    A0=qstab.chi(matrixParts(Q0));A1=qstab.chi(matrixParts(B));E0=eye(2*n);E1=zeros(2*n);
end
[L0,L1,L2,op]=qstab.internal.steinpencil(A0,A1,E0,E1);
R.diagnostics.operator=op;R.diagnostics.operatorDimension=op.dimension;
if any(~isfinite([L0(:);L1(:);L2(:)]))||rcond(L0)<eps*max(100,size(L0,1))
    R.message='The nominal generalized Stein operator is numerically singular or unrepresentable.';return
end
s0=norm(L0,2);s1=norm(L1,2);s2=norm(L2,2);
R.diagnostics.operatorCoefficientNorms=[s0,s1,s2];
logTau=Inf;
if s1>0,logTau=min(logTau,log(s0)-log(s1));end
if s2>0,logTau=min(logTau,.5*(log(s0)-log(s2)));end
if isinf(logTau)
    if isfinite(leadCandidate),R.message='Only a possibly degenerate leading contact remains unresolved.';return;end
    R=unbounded(R,'The generalized Stein operator is constant and nonsingular.');return
end
tau=exp(logTau);
if ~isfinite(tau)||tau==0,R.message='The internal parameter scale is unrepresentable.';return;end
Q0=L0/s0;Q1=zeros(size(L1));Q2=Q1;
if s1>0,Q1=(L1/s1)*exp(log(s1)+logTau-log(s0));end
if s2>0,Q2=(L2/s2)*exp(log(s2)+2*logTau-log(s0));end
if any(L1(:)~=0&Q1(:)==0)||any(L2(:)~=0&Q2(:)==0)
    R.message='Pencil scaling would erase a nonzero coefficient.';return
end
m=size(Q0,1);
if all(Q2(:)==0)
    PA=-Q0;PB=Q1;linearizationDegree=1;
else
    PA=[zeros(m),eye(m);-Q0,-Q1];PB=[eye(m),zeros(m);zeros(m),Q2];linearizationDegree=2;
end
try,[V,D,W]=builtin('eig',PA,PB,'qz');catch ME,R.message=['QZ failed: ',ME.message];return;end
values=diag(D);[eta,backward,pairing,bands]=eigenchecks(PA,PB,V,W,values);
finite=isfinite(values);realMask=finite&imag(values)==0;
possible=finite&(realMask|abs(imag(values))<=max(bands,10*eta));
reliable=realMask&isfinite(eta)&eta<=max(128*eps,o.Tolerance)*abs(values)&real(values)~=0;
kinds=repmat({'unitBoundary'},numel(values),1);
rawValues=values;rawEta=eta;
if isfinite(leadCandidate)
    leadScaled=leadCandidate/tau;
    if isfinite(leadScaled)&&leadScaled~=0
        values(end+1,1)=leadScaled;kinds{end+1,1}='degreeLoss';
        bands(end+1,1)=128*eps*abs(leadScaled);eta(end+1,1)=bands(end);
        possible(end+1,1)=true;reliable(end+1,1)=leadExact;
    end
end
[negative,negativeCandidate,negativeOK,ni]=endpoint(values,eta,bands,reliable,possible,-1);
[positive,positiveCandidate,positiveOK,pii]=endpoint(values,eta,bands,reliable,possible,1);
if any(isnan(values)),negativeOK=false;positiveOK=false;end
logTotalScale=logParameterScale+logTau;
R.candidateInterval=[scaleval(negativeCandidate,logTotalScale),scaleval(positiveCandidate,logTotalScale)];
R.diagnostics.scaledPencilEigenvalues=rawValues;
R.diagnostics.linearizationDegree=linearizationDegree;R.diagnostics.parameterScaleLog=logTotalScale;
R.diagnostics.eigenResidual=backward;R.diagnostics.eigenUncertaintyEstimate=rawEta;
R.diagnostics.finiteEigenpairPairing=pairing;R.diagnostics.infiniteEigenvalueCount=sum(isinf(rawValues));
R.diagnostics.leadingContact=leadCandidate;R.diagnostics.leadingContactExact=leadExact;
R.diagnostics.ambiguousScaledContacts=values(possible&~reliable);
R.diagnostics.boundaryVerification=NaN(2,3);
endpoints=[negative positive];okay=[negativeOK positiveOK];indices=[ni pii];
R.diagnostics.segmentGuards=cell(1,2);guardProcessed=0;guardTimer=tic;
for j=1:2
    if ~okay(j)&&isfinite(endpoints(j))&&indices(j)>0&&~any(isnan(values))
        targetBeta=endpoints(j)*tau*(1-o.Tolerance);
        S=segmentguard(L0,L1,L2,targetBeta,lead0,lead1, ...
            max(0,o.MaxGuardCells-guardProcessed),o.MaxGuardDepth,max(0,o.MaxGuardTime-toc(guardTimer)));
        R.diagnostics.segmentGuards{j}=S;guardProcessed=guardProcessed+S.processedCells;
        if S.accepted
            okay(j)=true;bounded=sort([scaleval(targetBeta,logParameterScale),scaleval(endpoints(j),logTotalScale)]);
            if j==1,R.negativeBracket=bounded;else,R.positiveBracket=bounded;end
        end
    end
    value=scaleval(endpoints(j),logTotalScale);event='none';witness=[];
    if okay(j)&&isfinite(endpoints(j))
        if ~isfinite(value)||value==0,okay(j)=false;
        else
            event=kinds{indices(j)};
            beta=endpoints(j)*tau;A=A0+beta*A1;E=E0+beta*E1;
            if strcmp(event,'unitBoundary')&&o.VerifyBoundary
                ev=builtin('eig',A,E,'qz');
                if any(~isfinite(ev)),okay(j)=false;
                else
                    R.diagnostics.boundaryVerification(j,:)=[min(abs(abs(ev)-1)),max(abs(ev)),rcond(E)];
                    okay(j)=min(abs(abs(ev)-1))<=max(1024*eps*m,o.Tolerance)&&max(abs(ev))<=1+max(1024*eps*m,o.Tolerance);
                end
            end
            if okay(j)
                [witness,okay(j)]=makeWitness(X0,X1,M,uv,value,side,typ,event,o.Tolerance);
            end
        end
    end
    if ~okay(j)
        if j==1,R.negativeBracket=[NaN NaN];else,R.positiveBracket=[NaN NaN];end
    end
    if okay(j)
        R.interval(j)=value;
        if j==1
            R.negative=value;R.negativeStatus='computed';R.negativeEvent=event;R.negativeWitness=witness;R.negativeParameterNorm=abs(value)*un;
        else
            R.positive=value;R.positiveStatus='computed';R.positiveEvent=event;R.positiveWitness=witness;R.positiveParameterNorm=abs(value)*un;
        end
    end
end
if all(okay),R.status='computed';R.message='Connected open Schur/regularity interval containing zero; numerical quadratic-pencil calculation.';
elseif any(okay),R.status='partial';R.message='Only one first-exit endpoint is resolved; earlier ambiguous contacts are not skipped.';
else,R.message='First contacts are numerically ambiguous; see candidateInterval.';end
end
function [A0,A1,E0,E1]=companionPair(C,B,n)
a0=zeros(n,n,4);a1=a0;e0=a0;e1=a0;
if n>1,a0(:,:,1)=diag(ones(n-1,1),1);e0(:,:,1)=diag([ones(1,n-1),0]);end
for k=1:n,a0(n,k,:)=reshape(-C(:,k),1,1,4);a1(n,k,:)=reshape(-B(:,k),1,1,4);end
e0(n,n,:)=reshape(C(:,end),1,1,4);e1(n,n,:)=reshape(B(:,end),1,1,4);
A0=qstab.chi(matrixParts(a0));A1=qstab.chi(matrixParts(a1));E0=qstab.chi(matrixParts(e0));E1=qstab.chi(matrixParts(e1));
end
function [value,exact,near]=leadingContact(a,b)
value=NaN;exact=false;near=false;
[bs,j]=max(abs(b));if bs==0,return;end
v=-a(j)/b(j);if ~isfinite(v)||v==0,return;end
r=a+v*b;scale=norm(a)+abs(v)*norm(b);
exact=all(r==0);near=norm(r)<=512*eps*max(scale,realmin);
if exact||near,value=v;end
end
function tf=fixedTriangular(A,B)
n=size(A,1);upper=true;lower=true;fixed=true;
for k=1:4
    aa=A(:,:,k);bb=B(:,:,k);
    upper=upper&&all(all(tril(aa,-1)==0))&&all(all(tril(bb,-1)==0));
    lower=lower&&all(all(triu(aa,1)==0))&&all(all(triu(bb,1)==0));
    fixed=fixed&&all(diag(bb)==0);
end
tf=n>0&&fixed&&(upper||lower);
end
function [W,ok]=makeWitness(X0,X1,M,u,alpha,side,typ,event,tol)
t=alpha*u;ok=all(isfinite(t));W=[];if ~ok,return;end
if strcmp(typ,'polynomial')
    Q=M.inputC0;D=zeros(size(Q));
    for k=1:M.grade+1
        if strcmp(side,'left'),D(:,k)=mul(t,M.inputC1(:,k));else,D(:,k)=mul(M.inputC1(:,k),t);end
        Q(:,k)=Q(:,k)+D(:,k);
    end
    ok=all(isfinite(Q(:)))&&all(isfinite(D(:)));if ~ok,return;end
    W=struct('kind',event,'valid',true,'alpha',alpha,'parameterComponents',t,'parameterNorm',norm(t), ...
        'coefficientComponents',Q,'perturbationComponents',D,'polynomial',polyParts(Q), ...
        'variable',M.variable,'testVariable','z','grade',M.grade,'side',side,'validated',false);
    if strcmp(event,'degreeLoss')
        idx=M.grade+1;if strcmp(M.variable,'zinv'),idx=1;end
        residual=norm(Q(:,idx))/max(realmin,norm(M.inputC0(:,idx))+norm(D(:,idx)));
        W.leadingResidual=residual;W.leadingIndexInInput=idx;ok=residual<=max(1024*eps,100*tol);W.valid=ok;
    end
else
    [a,b,c,d]=qstab.components(X0);Q=cat(3,a,b,c,d);
    [a,b,c,d]=qstab.components(X1);B=cat(3,a,b,c,d);D=zeros(size(B));
    for i=1:size(B,1),for j=1:size(B,2)
        q=reshape(B(i,j,:),4,1);if strcmp(side,'left'),q=mul(t,q);else,q=mul(q,t);end
        D(i,j,:)=reshape(q,1,1,4);
    end;end
    Q=Q+D;ok=all(isfinite(Q(:)));
    W=struct('kind',event,'valid',ok,'alpha',alpha,'parameterComponents',t,'parameterNorm',norm(t), ...
        'matrixComponents',Q,'matrix',matrixParts(Q),'side',side,'variable','z','validated',false);
end
end
function R=unbounded(R,message)
R.status='computed';R.interval=[-Inf Inf];R.candidateInterval=R.interval;
R.negative=-Inf;R.positive=Inf;R.negativeStatus='computed';R.positiveStatus='computed';
R.positiveParameterNorm=Inf;R.negativeParameterNorm=Inf;R.message=message;
end
function [eta,residual,pairing,band]=eigenchecks(A,B,V,W,values)
n=numel(values);eta=Inf(n,1);residual=NaN(n,1);pairing=NaN(n,1);band=NaN(n,1);
na=norm(A,2);nb=norm(B,2);arithmetic=64*eps*max(1,size(A,1));
for k=1:n
    if ~isfinite(values(k)),continue;end
    v=V(:,k);w=W(:,k);nv=norm(v);nw=norm(w);
    if ~(isfinite(nv)&&isfinite(nw)&&nv>0&&nw>0),continue;end
    v=v/nv;w=w/nw;sc=max(1,abs(values(k)));K=A/sc-(values(k)/sc)*B;
    back=max(norm(K*v),norm(K'*w));cs=na/sc+(abs(values(k))/sc)*nb;den=abs(w'*B*v);
    residual(k)=back/max(cs,realmin);pairing(k)=den;band(k)=arithmetic*sc;
    if den>0,eta(k)=(back+arithmetic*cs)/den*sc;end
end
end
function [value,candidate,ok,index]=endpoint(values,eta,band,reliable,possible,direction)
signed=direction*real(values);ids=find(reliable&signed>0);best=Inf;index=0;uncertainty=0;
if ~isempty(ids),[best,j]=min(signed(ids));index=ids(j);uncertainty=10*eta(index)+band(index);end
ids=find(possible&signed>0);if isempty(ids),candidate=direction*Inf;else,candidate=direction*min(signed(ids));end
ok=true;
for k=reshape(find(possible&~reliable),1,[])
    err=max(band(k),10*eta(k));
    if ~isfinite(err),ok=false;continue;end
    if signed(k)+err>0&&signed(k)-err<=best+uncertainty,ok=false;end
end
value=direction*best;
end
function y=scaleval(x,l)
if x==0,y=0;elseif isinf(x),y=x;elseif isnan(x),y=NaN;else,y=sign(x)*exp(log(abs(x))+l);end
end
function z=mul(x,y)
z=[x(1)*y(1)-x(2:4)'*y(2:4);x(1)*y(2:4)+y(1)*x(2:4)+cross(x(2:4),y(2:4))];
end
function P=polyParts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
function P=matrixParts(C)
P=struct('a',C(:,:,1),'b',C(:,:,2),'c',C(:,:,3),'d',C(:,:,4));
end

function S=segmentguard(L0,L1,L2,endpoint,lead0,lead1,maxCells,maxDepth,maxTime)
% A safe FIRST-contact guard: cover the interval from zero to just before
% a verified endpoint. Singular-value Lipschitz bounds rule out every earlier
% Stein singularity, including tangencies, without sampling acceptance.
S=struct('accepted',false,'processedCells',0,'acceptedCells',zeros(0,4), ...
    'parameterInterval',sort([0 endpoint]),'validated',false, ...
    'reason','Whole-segment nonsingularity not established.');
if maxCells<1||maxTime<=0||endpoint==0,return;end
lo=min(0,endpoint);hi=max(0,endpoint);timer=tic;m=size(L0,1);
if ~isempty(lead0)
    denom=lead1'*lead1;
    if denom==0,closest=0;else,closest=max(lo,min(hi,-(lead1'*lead0)/denom));end
    separation=norm(lead0+closest*lead1);
    if separation<=512*eps*(norm(lead0)+max(abs([lo hi]))*norm(lead1)),return;end
end
n0=norm(L0,2);n1=norm(L1,2);n2=norm(L2,2);queue=[lo hi 0];
while ~isempty(queue)
    if S.processedCells>=maxCells||toc(timer)>=maxTime,return;end
    node=queue(end,:);queue(end,:)=[];c=node(1)+(node(2)-node(1))/2;h=(node(2)-node(1))/2;
    H=L0+c*L1+c*c*L2;variation=h*norm(L1+2*c*L2,2)+h*h*n2;
    if any(~isfinite(H(:)))||~isfinite(variation),return;end
    sep=min(svd(H));guard=1024*eps*m*(n0+abs(c)*n1+c*c*n2+variation);
    S.processedCells=S.processedCells+1;
    if sep>variation+guard
        S.acceptedCells(end+1,:)=[node(1),node(2),node(3),sep-variation-guard]; %#ok<AGROW>
        continue
    end
    if node(3)>=maxDepth||c==node(1)||c==node(2),return;end
    queue(end+1,:)=[c,node(2),node(3)+1];queue(end+1,:)=[node(1),c,node(3)+1]; %#ok<AGROW>
end
S.accepted=true;S.reason='Full segment has a nonzero leading coefficient and nonsingular generalized Stein operator.';
end
function tf=integer(x,lo)
tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x==fix(x)&&x>=lo;
end
