function tests = qstab_cases_balls(writer)
%QSTAB_CASES_BALLS Weighted scalar coefficient balls, no private qfiles.
tests={ ...
 'Pointwise real first-order distances agree with all three Holder envelopes',@analyticReal; ...
 'A nonreal first-order centre minimizes over the full similarity class',@nonrealClass; ...
 'The algebraic class of s-j is not replaced by evaluation at i',@offSlice; ...
 'Exponent one allocates a critical perturbation to an optimal coefficient',@etaOne; ...
 'Exponent two reconstructs the Euclidean minimum-norm critical member',@etaTwo; ...
 'Exponent infinity uses independent full quaternion coefficient balls',@etaInfinity; ...
 'Noninteger Holder exponents preserve equality and witness residuals',@etaGeneral; ...
 'Zero-weight coefficients stay exactly fixed in returned polynomials',@fixedCoefficients; ...
 'At the zero class only the constant coefficient can cancel the value',@zeroArgument; ...
 'A zero uncertainty sensitivity gives infinity, not a fabricated witness',@zeroSensitivity; ...
 'Query array shape and conjugate class representatives are preserved',@queryShape; ...
 'Ascending vectors, tensors and scalar cells have the same ball semantics',@coefficientLayouts; ...
 'Common extreme coefficient and weight scalings preserve pointwise radii',@scaledInputs; ...
 'Homogeneous evaluation avoids overflow at large finite arguments',@largeArgument; ...
 'Tiny nonzero weights are not discarded as fixed coefficients',@smallWeights; ...
 'Malformed coefficient-ball weights are rejected',@invalidWeights; ...
 'Malformed Holder exponents are rejected',@invalidExponent; ...
 'Nonfinite queries and invalid pointwise options are rejected',@invalidArguments; ...
 'Underflow and unrepresentable critical polynomials remain explicit',@numericLimits; ...
 'Ball calculations do not alter inputs, directory or random state',@noSideEffects; ...
 'The CT first-order constant-coefficient margin has a global bracket',@ctFirst; ...
 'The full first-order l1, l2 and linfinity margins match closed forms',@ctFullNorms; ...
 'An infinite-frequency limiting event is leading-coefficient loss',@ctLeading; ...
 'The quadratic constant-coefficient margin is sqrt(3)/2',@ctQuadraticConstant; ...
 'The quadratic damping-coefficient ball has radius one',@ctQuadraticLinear; ...
 'A nonreal CT centre has its finite-frequency critical quaternion direction',@ctNoncommuting; ...
 'The Schur boundary includes both real endpoints',@schurReal; ...
 'The complete quaternion Schur sphere is covered, not just a fixed slice',@schurQuaternion; ...
 'Delay reversal reverses weights and maps the critical member back',@delayReverse; ...
 'Declared delay grade retains zero forward modes',@delayGrade; ...
 'A nonzero constant family loses the zero polynomial, not a finite mode',@constantBall; ...
 'All-zero weights leave the verified nominal family unchanged',@fixedBall; ...
 'Unstable, boundary and zero nominal members are distinguished',@nominalGuard; ...
 'A global cell budget preserves a nonconverged radius bracket',@cellBudgets; ...
 'Time and depth limits do not manufacture global convergence',@timeDepthBudgets; ...
 'Stored CT leaves cover the entire compactified frequency domain',@storedCoverage; ...
 'Bernstein matrix positivity is a whole-domain lower-radius certificate',@positiveCertificate; ...
 'Several analytic nonreal first-order margins stay inside computed brackets',@boundCrosscheck; ...
 'Common scaling preserves the global CT radius',@radiusScale; ...
 'Invalid variables, grades, tolerances and radius budgets are rejected',@radiusBadOptions; ...
 'Pointwise spectral-set decisions separate inside, outside and uncertainty',@setMembership; ...
 'A sampled spectral set makes no assertion about unsampled points',@setBoundaries; ...
 'A radius below the global lower bound gives strict robust stability',@stableDecision; ...
 'An admissible boundary-zero member disproves strict robust stability',@nonstableDecision; ...
 'Degree loss remains an assumption event in coefficient-ball decisions',@degreeDecision; ...
 'Unresolved radius brackets are never converted to instability',@unknownDecision; ...
 'Coefficient-ball data round-trip through the unchanged two-entry report ZIP',@()roundtrip(writer); ...
 'Selected qB02 ball examples return ordinary reproducible data',@exampleData};
end
function analyticReal()
for eta=[1 2 Inf]
    [r,I]=qstab.ballvalue([2 1],[1 .5],[0 1i 3i],'Exponent',eta);
    if eta==1,d=Inf;elseif isinf(eta),d=1;else,d=2;end
    exact=zeros(size(r));for j=1:numel(r),w=[0 1 3];exact(j)=sqrt(4+w(j)^2)/pn([1 .5*w(j)],d);end
    assert(all(I.valid)&&max(abs(r-exact))<1e-12);
end
end
function nonrealClass()
C=[2 1;0 0;3 0;0 0];z=[0 1 1i 3i -1+2i];
[r,I]=qstab.ballvalue(parts(C),[1 .5],z,'Witness',true);
exact=hypot(real(z)+2,abs(imag(z))-3)./sqrt(1+.25*abs(z).^2);
assert(all(I.valid)&&max(abs(r-exact))<1e-12);
for j=1:numel(z),checkWitness(parts(C),[1 .5],2,I.witnesses{j});end
end
function offSlice()
P=parts([0 1;0 0;-1 0;0 0]);[r,I]=qstab.ballvalue(P,[1 0],1i,'Witness',true);
assert(r<1e-14&&I.valid&&I.witnesses{1}.valid);
assert(norm(I.witnesses{1}.boundaryComponents-[0;0;1;0])<1e-12);
assert(norm(evaluate(qstab.polycomponents(P),[0;1;0;0]))>1);
end
function etaOne(),holderFixture(1);end
function etaTwo(),holderFixture(2);end
function etaInfinity(),holderFixture(Inf);end
function etaGeneral()
for eta=[1.1 1.5 3.5 1000],holderFixture(eta);end
end
function holderFixture(eta)
P=parts([2 .7 1; .1 .2 0; .4 -.3 0; -.2 .1 0]);g=[.3 .5 .2];
[r,I]=qstab.ballvalue(P,g,-.2+1.3i,'Exponent',eta,'Witness',true);
assert(I.valid&&I.lower<=r&&I.upper>=r);W=I.witnesses{1};checkWitness(P,g,eta,W);
assert(abs(W.weightedNorm-r)<1e-10*max(1,r));
if eta==1,assert(sum(any(W.perturbationComponents~=0,1))==1);end
end
function fixedCoefficients()
P=parts([2 3 1;.1 .2 0;0 -.3 0;0 0 0]);[~,I]=qstab.ballvalue(P,[1 0 0],.2+.8i,'Witness',true);
assert(I.valid);W=I.witnesses{1};C=qstab.polycomponents(P);
assert(isequal(W.coefficientComponents(:,2:3),C(:,2:3))&&W.fixedCoefficientsPreserved);
end
function zeroArgument()
[r,I]=qstab.ballvalue([2 3 1],[.5 2 3],0,'Witness',true);
assert(I.valid&&abs(r-4)<1e-12);W=I.witnesses{1};assert(all(W.coefficientComponents(:,1)==0));
assert(all(all(W.perturbationComponents(:,2:end)==0)));
end
function zeroSensitivity()
[r,I]=qstab.ballvalue([1 1],[0 1],0,'Witness',true);
assert(I.valid&&isinf(r)&&isinf(I.lower)&&isempty(I.witnesses{1}));
end
function queryShape()
z=reshape([1+2i 1-2i 0 .5],2,2);[r,I]=qstab.ballvalue([2 1],[1 0],z);
assert(isequal(size(r),size(z))&&r(1)==r(2)&&all(imag(I.classRepresentatives(:))>=0));
end
function coefficientLayouts()
P=[2+1i 3 1];T=reshape(P,1,1,3);C=num2cell(P);
a=qstab.ballvalue(P,[1 .2 0],1i);b=qstab.ballvalue(T,[1 .2 0],1i);c=qstab.ballvalue(C,[1 .2 0],1i);
assert(abs(a-b)<1e-13&&abs(a-c)<1e-13);
end
function scaledInputs()
r0=qstab.ballvalue([2 3 1],[1 .5 0],.3+1.2i);
for scale=[1e-200 1e200]
    [r,I]=qstab.ballvalue(scale*[2 3 1],scale*[1 .5 0],.3+1.2i,'Witness',true);
    assert(I.valid&&abs(r-r0)<1e-11*max(1,r0));
end
end
function largeArgument()
[r,I]=qstab.ballvalue([2 1],[1 1],1e150+2e150i);
assert(I.valid&&abs(r-1)<1e-12&&I.powerScaleLog>300);
end
function smallWeights()
[r,I]=qstab.ballvalue([1 1],[1e-250 0],0);
assert(I.valid&&abs(r/1e250-1)<1e-12);
end
function invalidWeights()
for g={[-1 0],[1],NaN,[1 Inf],{1,0},[1;0;1;0]}
    mustFail(@()qstab.ballvalue([1 1],g{1},1i));
end
end
function invalidExponent()
for e={0,.5,NaN,-Inf,[1 2],'2'}
    mustFail(@()qstab.ballvalue([1 1],[1 0],1i,'Exponent',e{1}));
end
end
function invalidArguments()
for z={[],NaN,Inf,1+Inf*1i,'abc'}
    mustFail(@()qstab.ballvalue([1 1],[1 0],z{1}));
end
mustFail(@()qstab.ballvalue([1 1],[1 0],1i,'Tolerance',0));
mustFail(@()qstab.ballvalue([1 1],[1 0],1i,'Witness',1));
end
function numericLimits()
[r,I]=qstab.ballvalue(1e-200*[1 1],1e200*[1 0],0,'Witness',true);
assert(~I.valid&&isnan(r)&&isempty(I.witnesses{1}));
[~,J]=qstab.ballvalue([2 1],[1 1],1e150+2e150i,'Witness',true);
assert(~J.valid||(~isempty(J.witnesses{1})&&J.witnesses{1}.evaluationResidual<1e-8));
mustFail(@()qstab.ballvalue([1e-300 1e300],[1 0],0));
end
function noSideEffects()
P=parts([2 1;.1 0;.3 0;-.2 0]);Q=P;g=[1 0];z=[0 1i];state=rng;wd=pwd;
qstab.ballvalue(P,g,z,'Witness',true);qstab.ballset(P,g,z,.2);
assert(isequal(P,Q)&&isequal(state,rng)&&strcmp(wd,pwd));
end
function ctFirst()
R=solve([2 1],[1 0]);bracket(R,2);assert(strcmp(R.limitingEvent,'finiteZero'));
checkWitness([2 1],[1 0],2,R.witness);
end
function ctFullNorms()
for eta=[1 2 Inf]
    R=solve([1 1],[1 1],'Exponent',eta);
    if isinf(eta),exact=1/sqrt(2);else,exact=1;end
    bracket(R,exact);
end
end
function ctLeading()
R=solve([1 1],[0 1]);bracket(R,1);
assert(strcmp(R.limitingEvent,'degreeLoss')&&isempty(R.witness.boundaryComponents));
assert(isequal(R.witness.coefficientComponents,[1 0;0 0;0 0;0 0]));
[~,S] = qstab.isstable(parts(R.witness.coefficientComponents(:,1)));assert(strcmp(S.status,'stable'));
end
function ctQuadraticConstant()
R=solve([1 1 1],[1 0 0]);bracket(R,sqrt(3)/2);
assert(strcmp(R.limitingEvent,'finiteZero'));checkWitness([1 1 1],[1 0 0],2,R.witness);
end
function ctQuadraticLinear()
R=solve([1 1 1],[0 1 0]);bracket(R,1);
end
function ctNoncommuting()
P=parts([2 1;0 0;3 0;0 0]);R=solve(P,[1 0]);bracket(R,2);
assert(abs(abs(R.boundaryArgument)-3)<.1);checkWitness(P,[1 0],2,R.witness);
end
function schurReal()
R=solve([-.5 1],[1 0],'Variable','z');bracket(R,.5);
T=solve([.5 1],[1 0],'Variable','z');bracket(T,.5);
assert(real(R.boundaryArgument)>0&&real(T.boundaryArgument)<0);
end
function schurQuaternion()
P=parts([-.3 1;0 0;-.4 0;0 0]);R=solve(P,[1 0],'Variable','z');bracket(R,.5);
checkWitness(P,[1 0],2,R.witness);
end
function delayReverse()
A=solve([-.5 1],[1 .2],'Variable','z');D=solve([1 -.5],[.2 1],'Variable','zinv');
assert(abs(A.lower-D.lower)<1e-12&&abs(A.upper-D.upper)<1e-12);
assert(isequal(D.witness.coefficientComponents,fliplr(A.witness.coefficientComponents)));
checkWitness([1 -.5],[.2 1],2,D.witness);
end
function delayGrade()
D=solve([1 -.5 0],[0 1 0],'Variable','zinv');bracket(D,.5);
assert(D.grade==2&&all(D.witness.coefficientComponents(:,3)==0));
checkWitness([1 -.5 0],[0 1 0],2,D.witness);
end
function constantBall()
[~,R] = qstab.ballradius(2+1i,2);bracket(R,sqrt(5)/2);
assert(R.gap==0&&strcmp(R.limitingEvent,'degreeLoss')&&all(R.witness.coefficientComponents(:)==0));
end
function fixedBall()
[~,R] = qstab.ballradius([1 1],[0 0]);assert(strcmp(R.status,'unbounded')&&isinf(R.lower)&&isinf(R.upper));
[~,S] = qstab.ballstable([1 1],[0 0],1e50);assert(strcmp(S.status,'stable'));
end
function nominalGuard()
[~,A] = qstab.ballradius([-1 1],[1 0]);[~,B] = qstab.ballradius([0 1],[1 0]);[~,Z] = qstab.ballradius([0 0],[1 0]);
assert(strcmp(A.status,'assumptionFailed')&&strcmp(B.status,'assumptionFailed')&&strcmp(Z.status,'assumptionFailed'));
assert(strcmp(A.nominal.status,'notStable')&&strcmp(B.nominal.status,'notStable')&&B.nominal.diagnostics.boundaryWitness.exactForSuppliedComponents);
end
function cellBudgets()
[~,R] = qstab.ballradius([1 1 1],[1 0 0],'AbsTol',1e-12,'RelTol',1e-12,'MaxCells',1,'MaxTime',Inf);
assert(strcmp(R.status,'bounded')&&R.diagnostics.processedCells==1);
assert(R.lower<=sqrt(3)/2&&R.upper>=sqrt(3)/2);
end
function timeDepthBudgets()
[~,R] = qstab.ballradius([1 1],[1 0],'MaxTime',0);assert(strcmp(R.status,'bounded')&&R.diagnostics.processedCells==0);
[~,D] = qstab.ballradius([-.5 1],[1 0],'Variable','z','MaxDepth',0,'MaxTime',Inf);
assert(strcmp(D.status,'bounded')&&D.lower<=.5&&D.upper>=.5);
end
function storedCoverage()
[~,R] = qstab.ballradius([1 1 1],[1 .2 0],'Exponent',3,'MaxCells',31,'MaxTime',Inf,'AbsTol',1e-12,'RelTol',1e-12);
T=sortrows(R.diagnostics.cells,1);assert(T(1,1)==0&&T(end,2)==1);
assert(abs(sum(T(:,2)-T(:,1))-1)<1e-14&&max(abs(T(1:end-1,2)-T(2:end,1)))<1e-14);
assert(R.diagnostics.processedCells<=31&&R.lower<=min(T(:,4))+1e-13);
end
function positiveCertificate()
R=solve([1 1],[1 1]);S=R.diagnostics.positivity;
assert(S.accepted&&S.attempted&&~S.validated&&all(S.acceptedCells(:,4)>0));
T=S.acceptedCells;assert(abs(sum(T(:,2)-T(:,1))-1)<1e-14);
end
function boundCrosscheck()
for imagPart=[.5 1 2]
    P=parts([1.5 1;0 0;imagPart 0;imagPart/2 0]);R=solve(P,[1 0]);bracket(R,1.5);
end
end
function radiusScale()
a=solve([1 1],[1 1]);
for scale=[1e-150 1e150]
    b=solve(scale*[1 1],scale*[1 1]);assert(abs(a.lower-b.lower)<1e-11&&abs(a.upper-b.upper)<1e-11);
end
end
function radiusBadOptions()
mustFail(@()qstab.ballradius([1 1],[1 0],'Variable','q'));
mustFail(@()qstab.ballradius([1 1],[1 0],'Grade',.5));
mustFail(@()qstab.ballradius([1 1],[1 0],'MaxCells',0));
mustFail(@()qstab.ballradius([1 1],[1 0],'MaxDepth',-1));
mustFail(@()qstab.ballradius([1 1],[1 0],'MaxTime',NaN));
mustFail(@()qstab.ballradius([1 1],[1 0],'AbsTol',0));
mustFail(@()qstab.ballradius([1 1],[1 0],'Tolerance',1));
end
function setMembership()
[~,S] = qstab.ballset([1 1],[1 0],[-1 -1+.1i 2],.2,'Witness',true);
assert(all(S.inside(1:2))&&S.outside(3)&&~any(S.unresolved));
end
function setBoundaries()
[~,S] = qstab.ballset([1 1],[1 0],0,1);
assert(S.unresolved&&S.isSampled&&~S.continuumCertified&&~S.validated);
end
function stableDecision()
[~,R] = qstab.ballstable([1 1],[1 0],.5,'AbsTol',1e-3,'RelTol',1e-3,'MaxTime',Inf);
assert(strcmp(R.status,'stable')&&R.isStable&&isempty(R.witness));
[~,Z] = qstab.ballstable([-.5 1],[1 0],.2,'Variable','z','AbsTol',1e-3,'RelTol',1e-3,'MaxTime',Inf);
assert(strcmp(Z.status,'stable'));
end
function nonstableDecision()
% The CLOSED unit ball contains s exactly at its boundary. One ulp
% inside, that witness is inadmissible even though it is in the band.
[tf,E]=qstab.ballstable([1 1],[1 0],1,'MaxTime',Inf);
assert(tf==0 && E.witness.exactBoundary.exactForSuppliedComponents);
[tf,E]=qstab.ballstable([1 1],[1 0],1-eps/2,'MaxTime',Inf);
assert(tf~=0 && ~strcmp(E.status,'notStable'));
[~,R] = qstab.ballstable([1 1],[1 0],1.2,'MaxTime',Inf);
assert(strcmp(R.status,'notStable')&&~R.isStable&&R.witness.valid&&R.witness.weightedNorm<=1.2);
[~,U] = qstab.ballstable([-1 1],[1 0],0);assert(strcmp(U.status,'notStable')&&strcmp(U.witness.kind,'nominal'));
end
function degreeDecision()
[~,R] = qstab.ballstable([1 1],[0 1],1.1,'MaxTime',Inf);
assert(strcmp(R.status,'assumptionFailed')&&isempty(R.isStable)&&strcmp(R.witness.kind,'degreeLoss'));
[~,D] = qstab.ballstable([1 0],[1 0],1.1,'Variable','zinv','MaxTime',Inf);
assert(strcmp(D.status,'assumptionFailed')&&all(D.witness.coefficientComponents(:,1)==0));
end
function unknownDecision()
[~,R] = qstab.ballstable([1 1 1],[1 0 0],.8,'MaxTime',0);
assert(strcmp(R.status,'undetermined')&&isempty(R.isStable)&&isempty(R.witness));
[~,Z] = qstab.ballstable([1 1],[1 0],0,'MaxTime',0);assert(strcmp(Z.status,'stable'));
mustFail(@()qstab.ballstable([1 1],[1 0],-1));
mustFail(@()qstab.ballset([1 1],[1 0],1i,Inf));
end
function roundtrip(writer)
root=tempname;mkdir(root);clean=onCleanup(@()rmdir(root,'s')); %#ok<NASGU>
E=run_qb02_balls('Verbose',false);
R=struct('toolbox',qstab.about(),'status','OVERALL_OK','ok',true,'passed',1,'failed',0,'skipped',0, ...
    'elapsedSeconds',0,'timestamp',datestr(now,30),'runId','ball-fixture', ...
    'environment',struct('version',version,'computer',computer,'quaternionPath','','qstabPath',''), ...
    'nativeMATLABRun',false,'nativeQuaternionAvailable',false,'nativeQuaternionOK',false, ...
    'cases',struct('group','fixture','name','Artificial coefficient-ball persistence fixture', ...
        'status','PASS','seconds',0,'message','','identifier',''), ...
    'coefficientBallExamples',E);
R=writer(R,root,fullfile(root,'raw'));dirout=fullfile(root,'readback');unzip(R.archiveFile,dirout);
A=load(fullfile(dirout,'report.mat'),'report');assert(isequaln(A.report.coefficientBallExamples,E));
items=dir(dirout);items=items(~[items.isdir]);assert(numel(items)==2);
end
function exampleData()
state=rng;wd=pwd;E=run_qb02_balls('Verbose',false);
assert(E.ok&&strcmp(E.schema,'qstab.ball-examples.1')&&isequal(state,rng)&&strcmp(wd,pwd));
assert(E.spectralSet.isSampled&&~E.spectralSet.continuumCertified);
end
function R=solve(P,g,varargin)
[~,R] = qstab.ballradius(P,g,'AbsTol',2e-3,'RelTol',2e-3,'MaxCells',2047,'MaxTime',Inf,varargin{:});
end
function bracket(R,exact)
assert(strcmp(R.status,'converged')&&~R.validated);
assert(R.lower<=exact+1e-10*max(1,exact)&&R.upper>=exact-1e-10*max(1,exact));
assert(R.gap<=.0021*(1+R.upper));
end
function checkWitness(P,g,eta,W)
assert(W.valid&&strcmp(W.kind,'finiteZero')&&W.evaluationResidual<1e-8);
C=qstab.polycomponents(P);Q=W.coefficientComponents;D=Q-C;
assert(isequal(D,W.perturbationComponents)&&all(all(D(:,g==0)==0)));
v=evaluate(Q,W.boundaryComponents);scale=sum(sqrt(sum(Q.^2,1)).*norm(W.boundaryComponents).^(0:size(Q,2)-1));
assert(norm(v)/max(realmin,scale)<1e-8);
a=zeros(size(g));for j=find(g>0),a(j)=norm(D(:,j))/g(j);end
assert(abs(pn(a,eta)-W.weightedNorm)<1e-10*max(1,W.weightedNorm));
end
function y=evaluate(C,x)
y=C(:,end);for j=size(C,2)-1:-1:1,y=mul(y,x)+C(:,j);end
end
function z=mul(x,y)
z=[x(1)*y(1)-x(2:4)'*y(2:4);x(1)*y(2:4)+y(1)*x(2:4)+cross(x(2:4),y(2:4))];
end
function P=parts(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
function y=pn(x,p)
x=abs(x);s=max(x);if s==0||isinf(p),y=s;else,y=s*sum((x/s).^p)^(1/p);end
end
function mustFail(f)
failed=false;try,f();catch,failed=true;end;assert(failed);
end
