function tests=qstab_cases_boxfix()
%QSTAB_CASES_BOXFIX Negative probes do not replace whole-family certificates.
tests={ ...
 'The original unstable CT box has an admissible witness with either backend',@originalCase; ...
 'The upper corner is checked when the lower corner alone is stable',@upperCorner; ...
 'Noncommuting fixed components survive a confirmed unstable box member',@noncommuting; ...
 'Forward and delay negative witnesses preserve the original declared grade',@forwardDelay; ...
 'A padded outside candidate is moved inward and its new spectrum is checked',@paddedCandidate; ...
 'Retraction that makes a candidate stable cannot preserve its old negative verdict',@stableRetraction; ...
 'Retraction onto the exact origin produces a strictly nonstable admissible member',@boundaryRetraction; ...
 'Outward-rounded original endpoints pass both strict component membership checks',@roundedEndpoint; ...
 'Fixed-component changes are rejected and candidate checks have no side effects',@fixedReject; ...
 'Stable box corners never certify a family with a hidden unstable edge member',@cornersNotCertificate; ...
 'A zero time budget performs no corner screening and cannot fabricate a witness',@zeroTime; ...
 'A sampled unstable member tightens only the upper radius bound at its actual norm',@upperOnly ...
};
end

function originalCase()
W=zeros(4,2);W(1,1)=1;
for method={'zeroexclusion','lyapunov'}
    [~,R] = qstab.boxstable([2 1],W,3,'Method',method{1},'MaxTime',Inf);
    checkWitness(R,[2 1;zeros(3,2)],W,3);
    assert(strcmp(R.witness.source,'originalLowerCorner'));
    assert(~R.diagnostics.backendCalled&&R.diagnostics.processedCells==0);
    assert(isequal(R.witness.coefficientComponents,[-1 1;zeros(3,2)]));
    assert(R.witness.spectralCheck.isConfirmedCounterexample);
end
end

function upperCorner()
W=zeros(4,2);W(1,1)=1;
[~,R] = qstab.boxstable([-2 -1],W,3,'MaxTime',Inf);
checkWitness(R,[-2 -1;zeros(3,2)],W,3);
assert(strcmp(R.witness.source,'originalUpperCorner'));
assert(numel(R.diagnostics.memberChecks)==2);
assert(strcmp(R.diagnostics.memberChecks{1}.status,'stable'));
end

function noncommuting()
C=[0 0;0 0;1 0;2 1];W=zeros(size(C));W(4,1)=1;
[~,R] = qstab.boxstable(parts(C),W,3,'MaxTime',Inf);
checkWitness(R,C,W,3);
assert(R.witness.coefficientComponents(3,1)==1);
assert(R.witness.coefficientComponents(4,2)==1);
assert(max(real(R.witness.spectralCheck.testRepresentatives))>.5);
end

function forwardDelay()
W=zeros(4,2);W(1,1)=1;
[~,R] = qstab.boxstable([-.5 1],W,1,'Variable','z','MaxTime',Inf);
checkWitness(R,[-.5 1;zeros(3,2)],W,1);
Wd=zeros(4,3);Wd(1,2)=1;
[~,D] = qstab.boxstable([1 -.5 0],Wd,1,'Variable','zinv','MaxTime',Inf);
checkWitness(D,[1 -.5 0;zeros(3,3)],Wd,1);
assert(D.grade==2&&all(D.witness.coefficientComponents(:,3)==0));
assert(any(D.witness.spectralCheck.testRepresentatives==0));
assert(any(isinf(D.witness.spectralCheck.inputRepresentatives)));
end

function paddedCandidate()
W=zeros(4,2);W(1,1)=1;M=qstab.internal.boxmodel([2 1],W,[]);
C=M.C;C(1,1)=-1-16*eps(5);
[V,D]=qstab.internal.boxcandidate(M,C,3,'s',1e-10);
assert(~isempty(V)&&V.valid&&V.insideRequestedBox);
assert(D.retracted&&V.weightedNorm<=3&&V.coefficientComponents(1,1)==-1);
assert(V.spectralCheck.isConfirmedCounterexample);
end

function stableRetraction()
W=zeros(4,2);W(1,1)=1;M=qstab.internal.boxmodel([2 1],W,[]);
C=M.C;C(1,1)=-1;
[V,D]=qstab.internal.boxcandidate(M,C,1,'s',1e-10);
assert(isempty(V)&&D.retracted&&D.insideRequestedBox);
assert(strcmp(D.status,'stable')&&D.coefficientComponents(1,1)==1);
end

function boundaryRetraction()
W=zeros(4,2);W(1,1)=1;M=qstab.internal.boxmodel([2 1],W,[]);
C=M.C;C(1,1)=-1;
[V,D]=qstab.internal.boxcandidate(M,C,2,'s',1e-10);
assert(~isempty(V)&&V.valid&&D.retracted&&D.insideRequestedBox);
assert(strcmp(D.status,'notStable')&&D.coefficientComponents(1,1)==0);
assert(D.spectralCheck.isConfirmedCounterexample && D.spectralCheck.nominal.diagnostics.boundaryWitness.exactForSuppliedComponents);
end

function roundedEndpoint()
W=zeros(4,2);W(1,1)=.2;M=qstab.internal.boxmodel([.1 1],W,[]);
C=M.C;C(1,1)=.1+.2;
assert(abs(C(1,1)-M.C(1,1))>W(1,1));
[V,D]=qstab.internal.boxcandidate(M,C,1,'s',1e-10);
assert(isempty(V)&&D.retracted&&D.inwardSteps>0&&D.insideRequestedBox);
assert(D.weightedNorm<=1&&abs(D.coefficientComponents(1,1)-.1)<=.2);
end

function fixedReject()
W=zeros(4,2);W(1,1)=1;M=qstab.internal.boxmodel([2 1],W,[]);
C=M.C;C(2,1)=1e-100;copyC=C;copyM=M;folder=pwd;state=rng;
[V,D]=qstab.internal.boxcandidate(M,C,3,'s',1e-10);
assert(isempty(V)&&~D.fixedComponentsPreserved&&~D.insideRequestedBox);
assert(isequal(C,copyC)&&isequal(M,copyM)&&strcmp(pwd,folder)&&isequal(rng,state));
end

function cornersNotCertificate()
% qB02: a square about (y,z)=(2.5,0), half-width 1. Its four corners
% satisfy y^2+z^2>3, but its edge midpoint (1.5,0) is unstable.
C=[1 1 1;2 0 0;0 2.5 0;0 0 0];W=zeros(size(C));W(3:4,2)=1;
[~,R] = qstab.boxstable(parts(C),W,1,'MaxCells',1,'MaxDepth',0,'MaxTime',Inf);
assert(strcmp(R.status,'undetermined')&&isempty(R.isStable));
assert(R.diagnostics.backendCalled&&numel(R.diagnostics.memberChecks)==2);
assert(all(cellfun(@(x)strcmp(x.status,'stable'),R.diagnostics.memberChecks)));
B=C;B(3,2)=1.5;B(4,2)=0;
[~,S] = qstab.isstable(parts(B),'Method','spectral');assert(strcmp(S.status,'notStable'));
end

function zeroTime()
W=zeros(4,2);W(1,1)=1;
[~,R] = qstab.boxstable([2 1],W,3,'MaxTime',0);
assert(strcmp(R.status,'undetermined')&&isempty(R.witness));
assert(isempty(R.diagnostics.memberChecks)&&R.diagnostics.processedCells==0);
end

function upperOnly()
C=[2 1;0 0;3 0;0 0];W=zeros(size(C));W(1,1)=1;
[~,R] = qstab.boxradius(parts(C),W,'InitialRadius',4,'MaxIterations',1,'MaxTime',Inf);
assert(R.lower==0&&isfinite(R.upper)&&R.upper>=2&&R.witness.valid);
assert(R.witness.insideRequestedBox&&R.witness.weightedNorm<=4);
assert(numel(R.diagnostics.trials)==1);
T=R.diagnostics.trials{1};assert(T.lowerBefore==T.lowerAfter&&T.upperAfter<T.upperBefore);
assert(~isempty(T.memberChecks)&&strcmp(T.status,'notStable'));
end

function checkWitness(R,C,W,r)
assert(strcmp(R.status,'notStable')&&isequal(R.isStable,false));
V=R.witness;assert(V.valid&&V.insideRequestedBox&&V.fixedComponentsPreserved);
delta=V.coefficientComponents-C;active=W>0;
assert(all(delta(~active)==0)&&all(abs(delta(active))<=r*W(active)));
assert(all(abs(delta(active))./W(active)<=r)&&V.weightedNorm<=r);
assert(V.spectralCheck.isConfirmedCounterexample);
end
function P=parts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
