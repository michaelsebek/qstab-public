function tests = qstab_cases_boxes(writer)
%QSTAB_CASES_BOXES Component directions, primal/dual gauges and full-box bounds.
tests={ ...
 'Rank-zero zonotopes distinguish zero targets from unattainable values',@rankZero; ...
 'A rank-one zonotope has the exact segment gauge and a primal representation',@rankOne; ...
 'A real interval has a finite rank-two gauge for every imaginary direction',@realPlane; ...
 'Rank-three and rank-four component gauges keep their intrinsic dimensions',@rankThreeFour; ...
 'A rotated rank-two zonotope does not project an out-of-span target into it',@rotatedRank; ...
 'Facet-budget exhaustion retains primal and dual bounds instead of infinity',@facetBudget; ...
 'Duplicated and zero generators do not change a zonotope gauge',@duplicates; ...
 'Box widths are component scales, not isotropic quaternion-ball radii',@notBalls; ...
 'An actual quaternion point is distinguished from its algebraic similarity class',@actualPoint; ...
 'Right and left evaluations reconstruct the respective noncommuting zeros',@sides; ...
 'Pure-imaginary compression uses at most eight generators at arbitrary degree',@eightGenerators; ...
 'Real query compression uses at most four component generators',@fourGenerators; ...
 'Zero-query cancellation is exact and leaves all higher coefficients unchanged',@zeroPoint; ...
 'Fixed real coefficient components remain bit-identical in witnesses',@fixedParts; ...
 'Nonzero tiny component widths are never silently declared fixed',@tinyWidth; ...
 'Common extreme coefficient and width scalings preserve pointwise box gauges',@scaling; ...
 'Vectors, tensors, coefficient cells and component width structs agree',@layouts; ...
 'Quaternion point arrays preserve their input shape',@queryShape; ...
 'Unrepresentable point evaluation returns unknown rather than a fabricated zero',@overflow; ...
 'Malformed widths and negative component scales are rejected',@invalidWidths; ...
 'Malformed evaluation options and facet budgets are rejected',@invalidPointOptions; ...
 'The qB02 hidden edge has a valid zonotope witness at its exact analytic radius',@edgePoint; ...
 'The qB02 full eight-dimensional example has its analytic boundary witness',@fullPoint; ...
 'A requested CT component box is decided as a whole family',@ctStable; ...
 'An unstable member is checked inside the original requested box',@ctUnstable; ...
 'Discrete and delay component-box decisions use the appropriate temporal variable',@dtStable; ...
 'Required leading-component loss is a separate assumption event',@leading; ...
 'Zero radius and all-zero widths retain the nominal family',@fixedFamily; ...
 'An unstable or ambiguous nominal input never gets a positive radius certificate',@badNominal; ...
 'The first-order real CT box has a global bracket around radius two',@ctRadius; ...
 'A real Schur box has a global bracket around radius one half',@dtRadius; ...
 'Delay reversal and declared zero padding preserve component-box radii',@delayRadius; ...
 'An explicit nonreal boundary seed gives an admissible upper bound, not a lower bound',@nonrealSeed; ...
 'A coefficient-box radius retains the original fixed components and witness norm',@radiusWitness; ...
 'No-iteration and zero-time budgets preserve honest one-sided or empty bounds',@zeroBudgets; ...
 'Cell budgets are global across all radius trials',@cellBudget; ...
 'A finite exploration cap is never presented as a finite critical radius',@radiusCap; ...
 'Unknown whole-family trials never tighten either side of a radius bracket',@unknownTrial; ...
 'Both established full-family methods can be used without renaming wrappers',@bothMethods; ...
 'Nonzero constant boxes lose the zero polynomial rather than a finite root',@constants; ...
 'Boundary seeds, time conventions and invalid radius options are checked',@invalidGlobalOptions; ...
 'Box calculations preserve inputs, current directory and random state',@sideEffects; ...
 'Component-box example payloads contain only ordinary serializable data',@exampleData; ...
 'Component-box results round-trip through the unchanged two-entry report ZIP',@()archive(writer) ...
};
end
function rankZero()
Z=qstab.internal.zonotopegauge(zeros(4,3),zeros(4,1),1e-10,100);assert(Z.valid&&Z.radius==0);
Z=qstab.internal.zonotopegauge(zeros(4,3),[1;0;0;0],1e-10,100);assert(Z.valid&&isinf(Z.radius)&&Z.infeasible);
end
function rankOne()
G=[1 2 -1;zeros(3,3)];x=[6;0;0;0];Z=qstab.internal.zonotopegauge(G,x,1e-10,100);
assert(Z.valid&&Z.numericalRank==1&&abs(Z.radius-1.5)<1e-10&&norm(G*Z.coefficients-x)<1e-9);
end
function realPlane()
W=zeros(4,3);W(1,:)=[.4 .6 0];
U=[1 0 0;0 1 0;0 0 1;1 2 3];
for k=1:size(U,1)
 u=U(k,:)'/norm(U(k,:));[r,I]=qstab.boxvalue([2 3 1],W,parts([0;2*u]),'Witness',true);
 assert(I.valid&&abs(r-5)<1e-8&&I.witnesses{1}.evaluationResidual<1e-8);
 assert(I.pointData{1}.gauge.numericalRank==2);
end
end
function rankThreeFour()
for d=3:4
 C=[1;2;3;4];C(d+1:4)=0;W=[.5;1;2;3];W(d+1:4)=0;
 [r,I]=qstab.boxvalue(parts(C),W,0,'Witness',true);
 assert(I.valid&&abs(r-2)<1e-10&&I.pointData{1}.gauge.numericalRank==d);
end
end
function rotatedRank()
G=[1 0;1 1;0 1;0 0];x=G*[2;-3];Z=qstab.internal.zonotopegauge(G,x,1e-10,100);
assert(Z.valid&&abs(Z.radius-3)<1e-9&&norm(G*Z.coefficients-x)<1e-8);
Z=qstab.internal.zonotopegauge(G,[0;0;0;1],1e-10,100);assert(Z.infeasible&&isinf(Z.radius));
Z=qstab.internal.zonotopegauge(G,[1;-1;1;0],1e-10,100);assert(~Z.valid||~isfinite(Z.upper));
end
function facetBudget()
G=[eye(4),eye(4)];x=[1;2;3;4];
Z=qstab.internal.zonotopegauge(G,x,1e-10,0);assert(Z.valid&&Z.facetCount==0);
assert(Z.lower<=2+1e-10&&Z.upper>=2-1e-10&&norm(G*Z.coefficients-x)<1e-8);
end
function duplicates()
G=[eye(4),zeros(4,1),-eye(4)];x=[2;4;6;8];
Z=qstab.internal.zonotopegauge(G,x,1e-10,4096);assert(Z.valid&&abs(Z.radius-4)<1e-9);
end
function notBalls()
P=parts([1;1;1;1]);[r,I]=qstab.boxvalue(P,ones(4,1),0,'Witness',true);
[b,B]=qstab.ballvalue(P,1,0,'Witness',true);assert(I.valid&&B.valid&&abs(r-1)<1e-10&&abs(b-2)<1e-10);
end
function actualPoint()
P=parts([0 1;0 0;-1 0;0 0]);W=zeros(4,2);W(1,1)=1;
[r,I]=qstab.boxvalue(P,W,parts([0;0;1;0]),'Witness',true);assert(I.valid&&r<1e-10);
[r,I]=qstab.boxvalue(P,W,1i);assert(isinf(r)&&I.valid&&~I.classMinimumClaim);
end
function sides()
C=[1 .3 1;.2 -.4 0;-.1 .2 0;.3 .1 0];W=ones(4,3);W(:,3)=0;q=[.2;-.4;.7;.1];
for side={'right','left'}
 [r,I]=qstab.boxvalue(parts(C),W,parts(q),'Evaluation',side{1},'Witness',true);
 assert(I.valid&&isfinite(r)&&I.witnesses{1}.valid&&I.witnesses{1}.evaluationResidual<1e-8);
 assert(norm(evalq(I.witnesses{1}.coefficientComponents,q,side{1}))<1e-8);
end
end
function eightGenerators()
C=zeros(4,13);C(1,:)=1:13;W=ones(size(C));W(:,end)=0;q=[0;1;2;3]/sqrt(14);
[r,I]=qstab.boxvalue(parts(C),W,parts(q),'Witness',true);
assert(I.valid&&isfinite(r)&&I.pointData{1}.generatorCount<=8&&I.witnesses{1}.valid);
end
function fourGenerators()
C=zeros(4,5);C(1,:)=[1 2 3 4 5];W=ones(size(C));W(:,end)=0;
[r,I]=qstab.boxvalue(parts(C),W,-2,'Witness',true);assert(I.valid&&isfinite(r)&&I.pointData{1}.generatorCount<=4);
end
function zeroPoint()
C=[2 1 .1;3 -.2 .4;4 .3 .2;5 -.4 1];W=ones(size(C));
[r,I]=qstab.boxvalue(parts(C),W,0,'Witness',true);Q=I.witnesses{1};
assert(I.valid&&abs(r-5)<1e-9&&all(Q.coefficientComponents(:,1)==0));
assert(isequal(Q.coefficientComponents(:,2:end),C(:,2:end))&&Q.evaluationResidual==0);
end
function fixedParts()
C=[2 1;0 .3;0 -.4;0 .5];W=zeros(size(C));W(1,1)=1;
[~,I]=qstab.boxvalue(parts(C),W,0,'Witness',true);Q=I.witnesses{1};
assert(I.valid&&all(Q.perturbationComponents(W==0)==0)&&Q.fixedComponentsPreserved);
end
function tinyWidth()
W=zeros(4,2);W(1,1)=1e-180;[r,I]=qstab.boxvalue([2 1],W,0,'Witness',true);
assert(I.valid&&isfinite(r)&&abs(r/2e180-1)<1e-11&&I.witnesses{1}.valid);
end
function scaling()
C=[2 1;.3 .1;-.4 0;.2 -.1];W=ones(size(C));W(:,end)=0;
[a,A]=qstab.boxvalue(parts(C),W,1i,'Witness',true);assert(A.valid);
for s=[1e-180 1e180]
 [b,B]=qstab.boxvalue(parts(s*C),s*W,1i,'Witness',true);assert(B.valid&&abs(a-b)<1e-8);
end
end
function layouts()
W=zeros(4,3);W(1,:)=[1 1 0];
[a,A]=qstab.boxvalue([1 1 1],W,1i,'Witness',true);
[b,B]=qstab.boxvalue(reshape([1 1 1],1,1,3),parts(W),1i,'Witness',true);
[c,C]=qstab.boxvalue({1,1,1},W,1i,'Witness',true);
assert(A.valid&&B.valid&&C.valid&&abs(a-b)<1e-12&&abs(a-c)<1e-12);
end
function queryShape()
W=ones(4,2);W(:,end)=0;z=reshape([0 1i 2 3i],2,2);[r,I]=qstab.boxvalue([2 1],W,z);
assert(isequal(size(r),size(z))&&isequal(size(I.valid),size(z))&&all(I.valid(:)));
end
function overflow()
W=ones(4,4);W(:,end)=0;[r,I]=qstab.boxvalue([1 1 1 1],W,1e300,'Witness',true);
assert(~I.valid&&isnan(r)&&I.lower==0&&isinf(I.upper));
end
function invalidWidths()
fails(@()qstab.boxvalue([1 1],[-1 0;zeros(3,2)],0));
fails(@()qstab.boxvalue([1 1],ones(1,2),0));
fails(@()qstab.boxvalue([1 1],complex(ones(4,2),ones(4,2)),0));
fails(@()qstab.boxvalue([1 1],NaN(4,2),0));
fails(@()qstab.boxvalue([1 1],ones(4,3),0));
end
function invalidPointOptions()
W=ones(4,2);fails(@()qstab.boxvalue([1 1],W,Inf));
fails(@()qstab.boxvalue([1 1],W,[]));fails(@()qstab.boxvalue([1 1],W,1i,'Evaluation','both'));
fails(@()qstab.boxvalue([1 1],W,1i,'MaxFacets',-.1));fails(@()qstab.boxvalue([1 1],W,1i,'Tolerance',0));
end
function edgePoint()
[C,W,q,ref]=edge();[r,I]=qstab.boxvalue(parts(C),W,parts(q),'Witness',true);
assert(I.valid&&abs(r-ref)<1e-8&&I.witnesses{1}.valid);
assert(all(I.witnesses{1}.perturbationComponents(W==0)==0));
end
function fullPoint()
[C,W,q,ref]=fullbox();[r,I]=qstab.boxvalue(parts(C),W,parts(q),'Witness',true);
assert(I.valid&&abs(r-ref)<1e-7&&I.witnesses{1}.evaluationResidual<1e-8);
end
function ctStable()
W=zeros(4,2);W(1,1)=1;[~,R] = qstab.boxstable([2 1],W,1);
assert(strcmp(R.status,'stable')&&R.isStable&&R.familyResult.certificateComplete);
end
function ctUnstable()
W=zeros(4,2);W(1,1)=1;[~,R] = qstab.boxstable([2 1],W,3);
assert(strcmp(R.status,'notStable')&&R.witness.valid&&R.witness.weightedNorm<=3);
assert(R.witness.fixedComponentsPreserved);
end
function dtStable()
W=zeros(4,2);W(1,1)=1;
[~,A] = qstab.boxstable([-.5 1],W,.2,'Variable','z');
[~,B] = qstab.boxstable([1 -.5],fliplr(W),.2,'Variable','zinv');
assert(strcmp(A.status,'stable')&&strcmp(B.status,'stable'));
end
function leading()
W=zeros(4,2);W(1,2)=1;[~,R] = qstab.boxstable([1 1],W,1);
assert(strcmp(R.status,'assumptionFailed')&&strcmp(R.witness.kind,'degreeLoss'));
assert(all(R.witness.coefficientComponents(:,2)==0));
W=zeros(4,2);W(1,1)=1;[~,R] = qstab.boxstable([1 -.5],W,1,'Variable','zinv');
assert(strcmp(R.status,'assumptionFailed')&&R.witness.leadingIndexInInput==1);
end
function fixedFamily()
[~,R] = qstab.boxstable([1 1],zeros(4,2),1e50);assert(strcmp(R.status,'stable'));
[~,R] = qstab.boxradius([1 1],zeros(4,2));assert(strcmp(R.status,'unbounded')&&isinf(R.lower)&&isinf(R.upper));
[~,R] = qstab.boxstable([1 1],ones(4,2),0);assert(strcmp(R.status,'stable'));
end
function badNominal()
W=ones(4,2);W(:,end)=0;
[~,R] = qstab.boxradius([-1 1],W);assert(strcmp(R.status,'assumptionFailed')&&~(R.lower>0));
[~,R] = qstab.boxstable([-1 1],W,0);assert(strcmp(R.status,'notStable'));
[~,R] = qstab.boxradius([0 1],W);assert(~(R.lower>0));
end
function ctRadius()
W=zeros(4,2);W(1,1)=1;[~,R] = qstab.boxradius([2 1],W,'RelTol',.02,'MaxCells',1024,'MaxTime',Inf);
bracket(R,2);assert(strcmp(R.status,'converged')&&R.diagnostics.globalLowerFromWholeFamilies);
end
function dtRadius()
W=zeros(4,2);W(1,1)=1;[~,R] = qstab.boxradius([-.5 1],W,'Variable','z','RelTol',.02,'MaxCells',2048,'MaxTime',Inf);
bracket(R,.5);assert(R.lower>0&&R.witness.valid);
end
function delayRadius()
W=zeros(4,3);W(1,2)=1;[~,R] = qstab.boxradius([1 -.5 0],W,'Variable','zinv','RelTol',.02,'MaxCells',2048,'MaxTime',Inf);
bracket(R,.5);assert(R.lower>0&&R.grade==2&&size(R.witness.coefficientComponents,2)==3);
end
function nonrealSeed()
C=[2 1;0 0;3 0;0 0];W=zeros(4,2);W(1,1)=1;
[~,R] = qstab.boxradius(parts(C),W,'BoundaryPoints',parts([0;0;-3;0]),'MaxIterations',0);
bracket(R,2);assert(R.lower==0&&R.witness.valid&&R.witness.evaluationResidual<1e-8);
assert(~R.diagnostics.directDirectionMinimumClaim);
end
function radiusWitness()
W=zeros(4,2);W(1,1)=1;[~,R] = qstab.boxradius([2 1],W,'MaxIterations',0);
Q=R.witness;assert(Q.valid&&Q.weightedNorm<=R.upper&&all(Q.perturbationComponents(W==0)==0));
end
function zeroBudgets()
W=zeros(4,2);W(1,1)=1;
[~,R] = qstab.boxradius([2 1],W,'MaxTime',0);assert(R.lower==0&&isinf(R.upper)&&R.diagnostics.iterations==0);
[~,R] = qstab.boxradius([2 1],W,'MaxIterations',0);bracket(R,2);assert(R.lower==0&&R.diagnostics.iterations==0);
end
function cellBudget()
W=zeros(4,2);W(1,1)=1;
[~,R] = qstab.boxradius([2 1],W,'MaxCells',3,'MaxTime',Inf);
assert(R.diagnostics.processedCells<=3&&R.diagnostics.iterations<=3);bracket(R,2);
end
function radiusCap()
W=zeros(4,2);W(1,1)=1;[~,R] = qstab.boxradius([2 1],W,'MaxRadius',.1,'MaxCells',512,'MaxTime',Inf);
assert(R.lower<=.1+1e-12&&R.upper>=2&&~strcmp(R.status,'converged'));
C=[2 1;0 0;3 0;0 0];
[~,R] = qstab.boxradius(parts(C),W,'MaxRadius',.1,'MaxCells',512,'MaxTime',Inf);
assert(R.lower<=.1+1e-12&&isinf(R.upper)&&~strcmp(R.status,'unbounded'));
end
function unknownTrial()
W=zeros(4,2);W(1,1)=1;
[~,R] = qstab.boxradius([-.5 1],W,'Variable','z','MaxDepth',0,'MaxIterations',4,'MaxTime',Inf);
for k=1:numel(R.diagnostics.trials)
 t=R.diagnostics.trials{k};
 if strcmp(t.status,'undetermined'),assert(t.lowerBefore==t.lowerAfter&&t.upperBefore==t.upperAfter);end
end
bracket(R,.5);
end
function bothMethods()
W=zeros(4,2);W(1,1)=1;
for method={'lyapunov','zeroexclusion'}
 [~,R] = qstab.boxradius([2 1],W,'Method',method{1},'RelTol',.05,'MaxCells',4096,'MaxTime',Inf);
 bracket(R,2);assert(R.lower>0);
end
end
function constants()
C=[1;2;3;4];W=[1;1;1;1];[~,R] = qstab.boxradius(parts(C),W);
assert(strcmp(R.status,'converged')&&R.lower==4&&R.upper==4&&strcmp(R.limitingEvent,'degreeLoss'));
[~,R] = qstab.boxstable(parts(C),W,3);assert(strcmp(R.status,'stable'));
[~,R] = qstab.boxstable(parts(C),W,4);assert(strcmp(R.status,'assumptionFailed'));
W=[1;0;0;0];[~,R] = qstab.boxradius(parts(C),W);assert(strcmp(R.status,'unbounded')&&isinf(R.lower));
end
function invalidGlobalOptions()
W=zeros(4,2);W(1,1)=1;
fails(@()qstab.boxradius([2 1],W,'BoundaryPoints',1));
fails(@()qstab.boxradius([-.5 1],W,'Variable','z','BoundaryPoints',2));
fails(@()qstab.boxradius([2 1],W,'Variable','other'));
fails(@()qstab.boxradius([2 1],W,'Method','sampling'));
fails(@()qstab.boxradius([2 1],W,'MaxCells',-1));
fails(@()qstab.boxradius([2 1],W,'InitialRadius',0));
fails(@()qstab.boxstable([2 1],W,-1));
fails(@()qstab.boxradius(ones(1,66),ones(4,66)));
fails(@()qstab.boxstable(ones(1,66),ones(4,66),.1));
end
function sideEffects()
P=[2 1];W=zeros(4,2);W(1,1)=1;sP=P;sW=W;d=pwd;s=rng;
qstab.boxvalue(P,W,1i);qstab.boxradius(P,W,'MaxIterations',0);qstab.boxstable(P,W,1);
assert(isequal(P,sP)&&isequal(W,sW)&&strcmp(pwd,d)&&isequal(rng,s));
end
function exampleData()
E=run_qb02_boxes('Verbose',false);assert(E.ok&&strcmp(E.schema,'qstab.component-box-examples.1'));
p=[tempname '.mat'];cleanup=onCleanup(@()delete(p)); %#ok<NASGU>
save(p,'E');T=load(p);assert(isequaln(E,T.E));
end
function archive(writer)
root=tempname;mkdir(root);cleanup=onCleanup(@()rmdir(root,'s')); %#ok<NASGU>
R=struct('toolbox',qstab.about(),'environment',struct('version',version,'computer',computer,'quaternionPath','','qstabPath',which('qstab.isstable')), ...
 'timestamp',datestr(now,30),'runId','box-fixture','nativeMATLABRun',false, ...
 'nativeQuaternionAvailable',false,'nativeQuaternionOK',false,'ok',true,'status','FIXTURE', ...
 'passed',1,'failed',0,'skipped',0,'elapsedSeconds',0, ...
 'cases',struct('group',{},'name',{},'status',{},'seconds',{},'message',{},'identifier',{}));
W=zeros(4,2);W(1,1)=1;[~,R.componentBoxExamples] = qstab.boxradius([2 1],W,'MaxIterations',0);
R=writer(R,root,fullfile(root,'raw'));dest=fullfile(root,'extract');mkdir(dest);names=unzip(R.archiveFile,dest);
assert(numel(names)==2);T=load(fullfile(dest,'report.mat'));assert(isequaln(T.report.componentBoxExamples,R.componentBoxExamples));
end
function bracket(R,x)
assert(R.lower<=x+1e-8&&R.upper>=x-1e-8&&~R.validated);
end
function [C,W,q,r]=edge()
C=[1 1 1;2 0 0;0 2.5 0;0 0 0];W=zeros(size(C));W(3:4,2)=1;
q=[0;-.5;0;-sqrt(3)/2];r=2.5-sqrt(3);
end
function [C,W,q,r]=fullbox()
C=[2 1 1;1 -1 0;-2 -1 0;1 1 0];W=ones(size(C));W(:,end)=0;r=(11-sqrt(103))/6;
a=(1-r)*[1;-1;-1;1];b=[2-r;1+r;-2-r;1+r];e=b-[4;0;0;0];o=2*a;
u=-mul([o(1);-o(2:4)],e)/(o'*o);q=[0;2*u(2:4)/norm(u(2:4))];
end
function v=evalq(C,q,side)
power=[1;0;0;0];v=zeros(4,1);
for k=1:size(C,2)
 if strcmp(side,'right'),v=v+mul(C(:,k),power);else,v=v+mul(power,C(:,k));end
 power=mul(power,q);
end
end
function z=mul(x,y),z=[x(1)*y(1)-x(2:4)'*y(2:4);x(1)*y(2:4)+y(1)*x(2:4)+cross(x(2:4),y(2:4))];end
function P=parts(C),P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]),'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));end
function fails(f),caught=false;try,f();catch,caught=true;end;assert(caught);end
