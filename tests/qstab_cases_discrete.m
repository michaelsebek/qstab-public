function tests = qstab_cases_discrete(writer)
%QSTAB_CASES_DISCRETE Analytic and adversarial ONE-parameter DT margins.
% Each case is new; the original 327 tests remain in their original files.
tests={ ...
 'Legacy CT radius subdivision keeps an immutable parent midpoint',@ctPartition; ...
 'Legacy CT radius witness frequency agrees with the returned boundary point',@ctWitnessFrequency; ...
 'Generalized Stein coefficients reproduce scalar norm squares',@steinScalar; ...
 'The generalized Stein operator acts on the full quaternion Hermitian space',@steinMatrix; ...
 'A scalar Schur ray includes both unit-circle endpoints',@rayFirst; ...
 'A quaternion state-update ray uses the discrete domain, not inverse eigenvalues',@rayMatrix; ...
 'An imaginary quaternion direction reaches a nonreal unit class',@rayImaginary; ...
 'The discrete direction is not normalized in the returned alpha interval',@rayDirection; ...
 'A tiny discrete direction retains distant finite contacts',@rayTiny; ...
 'Varying leading coefficients are kept affine without memberwise normalization',@rayLeading; ...
 'A proportional family reports zero-polynomial degree loss, not a unit root',@rayProportional; ...
 'Zero perturbations and fixed triangular diagonals give unbounded ray intervals',@rayInfinity; ...
 'A real Cayley transform preserves the first loss and later reentry',@rayReentry; ...
 'Noncommuting parameter sides keep their distinct first unit-circle contacts',@raySides; ...
 'An earlier defective unit-circle touch cannot be skipped for a later contact',@rayTouch; ...
 'Delay reversal preserves the ray parameter, grade and original witness coefficients',@rayDelay; ...
 'Common extreme scaling leaves discrete first-exit parameters unchanged',@rayScaling; ...
 'Discrete nominal failures and invalid ray options remain explicit',@rayInvalid; ...
 'Pointwise unit-class distances agree with a real first-order analytic formula',@pointReal; ...
 'The whole nonreal quaternion class is searched, not only a fixed complex slice',@pointNonreal; ...
 'Left and right pointwise witnesses cancel in the correct evaluation convention',@pointSides; ...
 'Unit endpoints use exact powers and preserve a vanishing perturbation',@pointEndpoints; ...
 'Discrete pointwise arrays, tensors and coefficient cells preserve shapes',@pointLayouts; ...
 'Common coefficient scaling preserves the unit-boundary quotient',@pointScaling; ...
 'Delay pointwise witnesses are reciprocated back to the input variable',@pointDelay; ...
 'Invalid angles, grades and unstable nominal inputs are rejected',@pointInvalid; ...
 'A constant perturbation has the exact first-order Schur global radius',@radiusFirst; ...
 'The parameter is common to all coefficients, not a coefficient ball',@radiusCorrelated; ...
 'Full leading-coefficient uncertainty retains the correct finite-boundary limit',@radiusLeading; ...
 'An identically zero polynomial is a degree-loss radius event',@radiusProportional; ...
 'A shifted nonreal unit-class minimum lies inside the global radius bracket',@radiusNonreal; ...
 'The full-circle Bernstein matrix shortcut proves a global lower bound',@radiusShortcut; ...
 'Unit-circle subdivision leaves a gap-free nonoverlapping partition',@radiusPartition; ...
 'Cell, depth and zero-time budgets retain conservative discrete radius brackets',@radiusBudgets; ...
 'Disabling stored angular cells does not disable the coverage check',@radiusNoStore; ...
 'Equivalent forward and delay families have equal one-parameter radius bounds',@radiusDelay; ...
 'Declared delay grade keeps a zero forward mode without losing its margin',@radiusGrade; ...
 'Zero perturbations and nonzero constants have explicit discrete radius semantics',@radiusConstants; ...
 'Non-Schur, boundary and degree-raising nominal cases never acquire a positive certificate',@radiusInvalidNominal; ...
 'One-parameter radii are invariant under a real Cayley change for both sides',@radiusCayley; ...
 'A real high-degree perturbation with endpoint zeros has a finite interior radius',@radiusInterior; ...
 'Invalid discrete radius options fail without altering the input family',@radiusInvalidOptions; ...
 'New discrete calculations preserve the working directory and random state',@noSideEffects; ...
 'Discrete namespace calls remain independent of private alias availability',@wrappers; ...
 'Discrete example payloads are ordinary reproducible data with checked witnesses',@examplePayload; ...
 'Discrete-margin data round-trip through the unchanged two-entry report ZIP',@()archivePayload(writer) ...
 };
end
function ctPartition()
P=parts([2 1;0 0;.37 0;0 0]);
for budget=[2 3 7]
 [~,R] = qstab.stabilityradius(P,1,'MaxCells',budget,'AbsTol',1e-12,'RelTol',1e-10);
 c=sortrows(R.diagnostics.cells,1);assert(c(1,1)==0&&c(end,2)==R.diagnostics.frequencyCutoff);
 assert(isequal(c(2:end,1),c(1:end-1,2))&&all(c(:,2)>c(:,1)));
 if budget==2,assert(c(1,2)==R.diagnostics.frequencyCutoff/2);end
end
end
function ctWitnessFrequency()
P=parts([2 1;0 0;.37 0;0 0]);
[~,R] = qstab.stabilityradius(P,1,'MaxCells',8,'AbsTol',1e-12,'RelTol',1e-10);
assert(norm(R.boundaryComponents(2:4))-R.frequency<1e-11&&abs(R.frequency-norm(R.boundaryComponents(2:4)))<1e-11);
[r,I]=qstab.freqmargin(P,1,R.frequency);assert(I.valid&&abs(r-R.upper)<1e-9);
end
function steinScalar()
a=[.2;.3;-.1;.4];b=[.1;-.2;.15;0];e=[1;.2;0;0];f=[0;.1;.1;.2];
[L0,L1,L2,I]=qstab.internal.steinpencil(Q(a),Q(b),Q(e),Q(f));
close(L0,norm(e)^2-norm(a)^2,1e-13);close(L1,2*(e'*f-a'*b),1e-13);close(L2,norm(f)^2-norm(b)^2,1e-13);
assert(I.dimension==1);
end
function steinMatrix()
A=[.1 .2;0 -.3];B=[0 .4;0 .2];E=eye(2);F=[.05 0;0 -.1];
[L0,L1,L2,I]=qstab.internal.steinpencil(qstab.chi(A),qstab.chi(B),qstab.chi(E),qstab.chi(F));
assert(I.dimension==6&&isreal(L0)&&isreal(L1)&&isreal(L2));
% X=I has coordinates [1;1;0;0;0;0]. A real symmetric image has no i,j,k parts.
for t=[-.4 0 .7]
    Y=(E+t*F)*(E+t*F)'-(A+t*B)*(A+t*B)';expected=[diag(Y);sqrt(2)*Y(1,2);0;0;0];
    assert(norm((L0+t*L1+t*t*L2)*[1;1;0;0;0;0]-expected)<1e-12);
end
end
function rayFirst()
[~,R] = qstab.raymargin([-.5 1],1,1,'Variable','z');computed(R);close(R.negative,-.5,1e-9);close(R.positive,1.5,1e-9);
assert(strcmp(R.positiveEvent,'unitBoundary'));member(R.positiveWitness);
end
function rayMatrix()
[~,R] = qstab.raymargin(.5,1,1,'Type','matrix','Domain','discrete');computed(R);close(R.negative,-1.5,1e-9);close(R.positive,.5,1e-9);
assert(strcmp(R.variable,'z')&&norm(R.positiveWitness.matrix.a-1)<1e-9);
end
function rayImaginary()
[~,R] = qstab.raymargin([-.5 1],1,parts([0;1;0;0]),'Variable','z');computed(R);close(R.positive,sqrt(3)/2,1e-9);close(R.negative,-sqrt(3)/2,1e-9);
end
function rayDirection()
[~,R] = qstab.raymargin([-.5 1],1,2,'Variable','z');computed(R);close(R.negative,-.25,1e-9);close(R.positive,.75,1e-9);close(R.positiveParameterNorm,1.5,1e-9);
end
function rayTiny()
[~,R] = qstab.raymargin([-.5 1],1,1e-16,'Variable','z');computed(R);close(R.negative,-.5e16,1e-9);close(R.positive,1.5e16,1e-9);
end
function rayLeading()
[~,R] = qstab.raymargin([-.5 1],[0 1],1,'Variable','z');computed(R);close(R.negative,-.5,1e-9);assert(isinf(R.positive)&&R.positive>0);
assert(strcmp(R.negativeEvent,'unitBoundary'));
end
function rayProportional()
[~,R] = qstab.raymargin([-.5 1],[-.5 1],1,'Variable','z');computed(R);close(R.negative,-1,1e-12);assert(isinf(R.positive));assert(strcmp(R.negativeEvent,'degreeLoss'));
assert(R.negativeWitness.leadingResidual<1e-12);
end
function rayInfinity()
[~,R] = qstab.raymargin([-.5 1],0,1,'Variable','z');computed(R);assert(all(isinf(R.interval)));
[~,R] = qstab.raymargin([.2 .7;0 -.4],[0 3;0 0],1,'Type','matrix','Domain','discrete');computed(R);assert(all(isinf(R.interval)));
end
function rayReentry()
P=cayley([13/16 1 1 1],3);B=cayley([3 1 1],3);
[~,R] = qstab.raymargin(P,B,1,'Variable','z');computed(R);close(R.positive,1/4,2e-8);close(R.negative,-13/48,2e-8);
end
function raySides()
[P,B,u]=sideData();
for right=[false true]
    side='left';reference=1;if right,side='right';reference=2-sqrt(3);end
    [~,R] = qstab.raymargin(cayley(P,2),cayley(B,2),u,'Side',side,'Variable','z');
    assert(strcmp(R.positiveStatus,'computed'));close(R.positive,reference,2e-8);member(R.positiveWitness);
    if strcmp(R.negativeStatus,'computed'),assert(isinf(R.negative)&&R.negative<0);else,assert(isnan(R.negative));end
end
end
function rayTouch()
P=parts([1 2 1;zeros(1,3);1 0 0;zeros(1,3)]);B=parts([-1 -1;zeros(1,2);-1 -1;zeros(1,2)]);
[~,R] = qstab.raymargin(cayley(P,2),cayley(B,2),1,'Variable','z');
if strcmp(R.positiveStatus,'computed'),close(R.positive,1,2e-6);
else,assert(isnan(R.positive)&&any(strcmp(R.status,{'partial','undetermined'})));end
assert(~(strcmp(R.positiveStatus,'computed')&&abs(R.positive-1.5)<1e-4));
end
function rayDelay()
[~,R] = qstab.raymargin([1 -.5 0],[0 1 0],1,'Variable','zinv');computed(R);close(R.negative,-.5,1e-9);close(R.positive,1.5,1e-9);
assert(R.grade==2&&strcmp(R.positiveWitness.variable,'zinv'));member(R.positiveWitness);
end
function rayScaling()
for sc=[1e-180,1e180,-3]
    [~,R] = qstab.raymargin(sc*[-.5 1],sc,1,'Variable','z');computed(R);close(R.negative,-.5,1e-9);close(R.positive,1.5,1e-9);
end
end
function rayInvalid()
[~,R] = qstab.raymargin([-2 1],1,1,'Variable','z');assert(strcmp(R.status,'assumptionFailed'));
mustFail(@()qstab.raymargin([-.5 1],1,0,'Variable','z'));
mustFail(@()qstab.raymargin([-.5 1],1,1,'Variable','z','Grade',0));
mustFail(@()qstab.raymargin(.5,1,1,'Type','matrix','Domain','discrete','Variable','zinv'));
mustFail(@()qstab.raymargin([-.5 1],1,1,'Variable','z','MaxGuardCells',-1));
end
function pointReal()
t=[0 pi/3 pi/2 pi];[r,I]=qstab.freqmargin([-.5 1],1,t,'Variable','z');
assert(all(I.valid));assert(max(abs(r-sqrt(1.25-cos(t))))<1e-10);
end
function pointNonreal()
a=[.2;0;.3;0];P=parts([-a,[1;0;0;0]]);theta=atan2(.3,.2);
[r,I]=qstab.freqmargin(P,1,theta,'Variable','z');assert(I.valid);close(r,1-norm(a),1e-10);
assert(abs(I.boundaryComponents(3))>.5&&abs(I.boundaryComponents(2))<1e-8);
end
function pointSides()
[P,B,~]=sideData();P=cayley(P,2);B=cayley(B,2);
for side={'left','right'}
    [~,I]=qstab.freqmargin(P,B,[.3 .9 1.7],'Side',side{1},'Variable','z');assert(all(I.valid));
    for k=1:3,assert(I.witnesses{k}.evaluationResidual<1e-8);member(I.witnesses{k});end
end
end
function pointEndpoints()
[r,I]=qstab.freqmargin([0 0 1],[-1 0 1],[0 pi],'Variable','z');assert(all(isinf(r))&&all(I.valid));
[r,I]=qstab.freqmargin([.5 1],1,[0 pi],'Variable','z');assert(all(I.valid));close(r(1),1.5,1e-10);close(r(2),.5,1e-10);
end
function pointLayouts()
[r,I]=qstab.freqmargin(reshape([-.5 1],1,1,2),{1},[0;1],'Variable','z');assert(isequal(size(r),[2 1])&&isequal(size(I.valid),[2 1]));
q=qstab.freqmargin({-.5,1},1,[0;1],'Variable','z');assert(norm(r-q)<1e-12);
end
function pointScaling()
base=qstab.freqmargin([-.5 1],1,[.2 1.2],'Variable','z');
for sc=[1e-180 1e180]
    r=qstab.freqmargin(sc*[-.5 1],sc,[.2 1.2],'Variable','z');assert(norm(r-base)<1e-10);
end
end
function pointDelay()
P=parts([-.2 1;0 0;-.3 0;0 0]);Pd=reverse(P);B=1;Bd=[0 1];
[r,I]=qstab.freqmargin(P,B,.8,'Variable','z');[s,J]=qstab.freqmargin(Pd,Bd,.8,'Variable','zinv');close(r,s,1e-12);
assert(norm(I.parameterComponents-J.parameterComponents)<1e-10);
assert(norm(I.boundaryComponents-[J.boundaryComponents(1);-J.boundaryComponents(2:4)])<1e-10);member(J.witnesses{1});
end
function pointInvalid()
mustFail(@()qstab.freqmargin([-.5 1],1,pi+.1,'Variable','z'));
mustFail(@()qstab.freqmargin([-.5 1],1,-1,'Variable','z'));
mustFail(@()qstab.freqmargin([-2 1],1,0,'Variable','z'));
mustFail(@()qstab.freqmargin([-.5 1],[0 0 1],0,'Variable','z'));
end
function radiusFirst()
R=margin([-.5 1],1);bracket(R,.5);assert(strcmp(R.limitingEvent,'finiteZero'));member(R.witness);
end
function radiusCorrelated()
R=margin([-.25 1],[.5 1]);bracket(R,.5);assert(R.diagnostics.completeDomainCovered);
end
function radiusLeading()
R=margin([-.5 1],[0 1]);bracket(R,.5);close(R.degreeLossRadius,1,1e-12);assert(strcmp(R.limitingEvent,'finiteZero'));
end
function radiusProportional()
R=margin([-.5 1],[-.5 1]);bracket(R,1);assert(strcmp(R.limitingEvent,'degreeLoss'));assert(R.witness.leadingResidual<1e-9);
end
function radiusNonreal()
a=[.2;0;.3;0];R=margin(parts([-a,[1;0;0;0]]),1);bracket(R,1-norm(a));member(R.witness);
end
function radiusShortcut()
R=margin([-.5 1],1);assert(R.diagnostics.bernsteinShortcut.accepted);assert(R.diagnostics.completeDomainCovered);
end
function radiusPartition()
P=parts([-.2 1;0 0;-.3 0;0 0]);
for budget=[1 3 7 15]
 [~,R] = qstab.stabilityradius(P,1,'Variable','z','MaxCells',budget,'RelTol',1e-8,'AbsTol',1e-10);
 assert(R.lower<=1-sqrt(.13)+1e-9&&R.upper>=1-sqrt(.13)-1e-9);assert(R.diagnostics.partitionVerified);
 if ~R.diagnostics.bernsteinShortcut.accepted
  c=R.diagnostics.cells;assert(c(1,1)==0&&c(end,2)==pi&&isequal(c(2:end,1),c(1:end-1,2)));
 end
end
end
function radiusBudgets()
P=parts([-.2 1;0 0;-.3 0;0 0]);
for args={{'MaxCells',1},{'MaxCells',1,'MaxDepth',0},{'MaxTime',0}}
 argsNow=args{1};[~,R] = qstab.stabilityradius(P,1,'Variable','z',argsNow{:},'AbsTol',1e-12,'RelTol',1e-10);
 assert(any(strcmp(R.status,{'bounded','undetermined'}))&&R.lower<=1-sqrt(.13)+1e-9);assert(R.diagnostics.partitionVerified);
end
end
function radiusNoStore()
P=parts([-.2 1;0 0;-.3 0;0 0]);
[~,A] = qstab.stabilityradius(P,1,'Variable','z','MaxCells',7,'StoreCells',true);
[~,B] = qstab.stabilityradius(P,1,'Variable','z','MaxCells',7,'StoreCells',false);
assert(A.lower==B.lower&&A.upper==B.upper&&isempty(B.diagnostics.cells)&&B.diagnostics.partitionVerified);
end
function radiusDelay()
A=margin([-.5 1],1);B=margin([1 -.5],[0 1],'Variable','zinv');bracket(B,.5);
assert(abs(A.lower-B.lower)<1e-10&&abs(A.upper-B.upper)<1e-10);member(B.witness);
end
function radiusGrade()
R=margin([1 -.5],[0 1],'Variable','zinv','Grade',2);bracket(R,.5);assert(R.grade==2&&size(R.witness.coefficientComponents,2)==3);
end
function radiusConstants()
R=margin([-.5 1],0);assert(strcmp(R.status,'unbounded')&&isinf(R.lower));
R=margin(2,1);bracket(R,2);assert(strcmp(R.limitingEvent,'degreeLoss'));
end
function radiusInvalidNominal()
for P={[-2 1],[1 1],[0 0]}
 R=margin(P{1},1);assert(any(strcmp(R.status,{'assumptionFailed','undetermined'})));
end
R=margin([-.5 1],[0 0 1]);assert(strcmp(R.status,'assumptionFailed'));
end
function radiusCayley()
[P,B,~]=sideData();
for side={'left','right'}
 [~,A] = qstab.stabilityradius(P,B,'Side',side{1},'AbsTol',1e-5,'RelTol',.01,'MaxCells',1023);
 R=margin(cayley(P,2),cayley(B,2),'Side',side{1},'MaxCells',4095);
 assert(isfinite(A.upper)&&isfinite(R.upper));assert(A.lower<=R.upper+1e-8&&R.lower<=A.upper+1e-8);
 assert(strcmp(R.status,'converged'));member(R.witness);
end
end
function radiusInterior()
R=margin([0 0 0 0 1],[-1 0 0 0 1],'MaxCells',4095);bracket(R,.5);
end
function radiusInvalidOptions()
mustFail(@()margin([-.5 1],1,'Variable','q'));
mustFail(@()margin([-.5 1],1,'MaxCells',0));
mustFail(@()margin([-.5 1],1,'MaxDepth',-1));
mustFail(@()margin([-.5 1],1,'Tolerance',1));
mustFail(@()margin([-.5 1],1,'Grade',0));
end
function noSideEffects()
pwd0=pwd;state=rng;P=parts([-.2 1;0 0;-.3 0;0 0]);saved=P;
qstab.freqmargin(P,1,.7,'Variable','z');margin(P,1);qstab.raymargin(P,1,1,'Variable','z');
assert(strcmp(pwd,pwd0)&&isequal(rng,state)&&isequal(P,saved));
end
function wrappers()
[~,A] = qstab.raymargin([-.5 1],1,1,'Variable','z');
assert(max(abs(A.interval-[-.5 1.5]))<1e-8);
[~,B] = qstab.stabilityradius([-.5 1],0,'Variable','z');
assert(strcmp(B.status,'unbounded'));
end
function examplePayload()
D=run_qb01_discrete('Verbose',false);assert(D.ok&&isstruct(D));
end
function archivePayload(writer)
root=tempname;mkdir(root);cleanup=onCleanup(@()rmdir(root,'s')); %#ok<NASGU>
R=struct('toolbox',qstab.about(),'environment',struct('version',version,'computer',computer,'quaternionPath','','qstabPath',which('qstab.isstable')), ...
 'runId','discrete-roundtrip','timestamp','unit','nativeMATLABRun',false,'nativeQuaternionOK',false, ...
 'status','OVERALL_OK','passed',1,'failed',0,'skipped',0,'elapsedSeconds',0, ...
 'cases',struct('status','PASS','group','discrete','name','archive','seconds',0,'message',''));
[~,ray]=qstab.raymargin([-.5 1],1,1,'Variable','z');
R.discreteMarginExamples=struct('radius',margin([-.5 1],1),'ray',ray);
R=writer(R,root,fullfile(root,'qstab_test_results'));folder=fullfile(root,'unzipped');unzip(R.archiveFile,folder);S=load(fullfile(folder,'report.mat'),'report');
assert(S.report.discreteMarginExamples.radius.lower==R.discreteMarginExamples.radius.lower&&strcmp(S.report.discreteMarginExamples.ray.status,'computed'));
end
function R=margin(P,B,varargin)
[~,R] = qstab.stabilityradius(P,B,'Variable','z','AbsTol',1e-5,'RelTol',.005,'MaxCells',2047,varargin{:});
end
function bracket(R,x)
assert(strcmp(R.status,'converged')&&R.lower<=x+1e-8*max(1,x)&&R.upper>=x-1e-8*max(1,x));
assert(R.upper-R.lower<=1e-5+.005*R.upper+1e-10);
end
function member(W)
assert(W.valid);
if isfield(W,'evaluationResidual'),assert(W.evaluationResidual<1e-8);end
if isfield(W,'coefficientComponents')&&strcmp(W.kind,'finiteZero')
 [z,I]=qstab.roots(W.polynomial,'Grade',W.grade);assert(I.isValid&&min(abs(abs(z(isfinite(z)))-1))<1e-5);
elseif isfield(W,'coefficientComponents')&&strcmp(W.kind,'unitBoundary')
 [z,I]=qstab.roots(W.polynomial,'Grade',W.grade);assert(I.isValid&&min(abs(abs(z(isfinite(z)))-1))<1e-5);
end
end
function computed(R),assert(strcmp(R.status,'computed')&&strcmp(R.positiveStatus,'computed')&&strcmp(R.negativeStatus,'computed'));end
function close(a,b,t),assert(isfinite(a)&&abs(a-b)<=t*max(1,abs(b)));end
function mustFail(f)
ok=false;try,f();catch,ok=true;end;assert(ok,'Expected a rejected input.');
end
function P=parts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]),'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
function B=Q(v)
B=[v(1)+1i*v(2),v(3)+1i*v(4);-v(3)+1i*v(4),v(1)-1i*v(2)];
end
function P=reverse(P),C=qstab.polycomponents(P);P=parts(fliplr(C));end
function P=cayley(P,n)
C=qstab.polycomponents(P);R=zeros(4,n+1);
for k=0:size(C,2)-1
 b=1;for j=1:k,b=conv(b,[-1 1]);end;for j=1:n-k,b=conv(b,[1 1]);end
 R=R+C(:,k+1)*b;
end
P=parts(R);
end
function [P,B,u]=sideData()
P=parts([1 2 1;0 0 0;1 0 0;0 0 0]);B=parts([0 0;0 1;0 0;1 0]);u=parts([0;1;0;1]);
end
