function cases = qstab_cases_diagnostics(writer)
%QSTAB_CASES_DIAGNOSTICS Read-only geometry, provenance and witness tests.
cases={ ...
 'Bernstein controls include mixed terms, not only vertex norms',@bernstein; ...
 'A coarse outer hull containing zero is not an instability witness',@falseAlarm; ...
 'Subdivision keeps separate local hulls and original convex weights',@splitHull; ...
 'Explicit subcells do not silently claim full-domain coverage',@coverageFlag; ...
 'Real, point and collinear hulls are retained without Qhull assumptions',@degenerateHull; ...
 'Independent imaginary-axis boxes have the correct norm-coordinate ranges',@evenOdd; ...
 'Shared affine coefficients are not independently recombined',@correlated; ...
 'Discrete arguments use a joint quadratic enclosure',@unitCircle; ...
 'Value enclosures in delay notation evaluate the original input variable',@inputVariable; ...
 'Sub-boxes preserve the exact original parameter indexing',@parameterBoxes; ...
 'Noncommuting parameter sides give different reconstructed images',@sides; ...
 'Invalid weights, boxes and nonfinite arguments are rejected',@invalidInputs; ...
 'Enclosure overflow is explicit, not a fabricated zero or empty set',@overflow; ...
 'Diagnostics do not alter inputs, directory or random state',@noSideEffects; ...
 'Retained Bernstein cells keep the reported solver status and storage count',@retained; ...
 'A stable result without stored leaves remains stable but storage-incomplete',@missingStorage; ...
 'Unvisited subfamilies remain explicit non-frequency obligations',@pending; ...
 'Lyapunov slack is labelled separately from norm-value separation',@lyapunov; ...
 'Constant-family completion does not manufacture stored leaves',@constant; ...
 'A saved Bernstein cell reconstructs a whole-frequency outer enclosure',@certificatePoly; ...
 'Four-generator certificate weights map back to the complete input list',@certificateFour; ...
 'A mismatched family cannot silently reconstruct a saved Bernstein cell',@wrongFamily; ...
 'Saved interval and affine rectangles need no robust solver rerun',@certificateRectangles; ...
 'Delay certificate geometry explicitly uses forward z',@certificateDelay; ...
 'Frequency disk rendering is circumscribed rather than inscribed',@outerDisk; ...
 'Lyapunov and malformed certificates have no invented norm-value polygon',@noGeometry; ...
 'An admissible unstable polytope member passes an independent spectral check',@witnessPoly; ...
 'Tampered coefficients or convex weights cannot pass membership checks',@witnessTampered; ...
 'Interval witnesses are checked against every original component bound',@witnessInterval; ...
 'Affine witness reconstruction preserves shared real parameters',@witnessAffine; ...
 'Delay witnesses retain both reciprocal input and forward representatives',@witnessDelay; ...
 'A zero forward mode is retained as an infinite delay representative',@witnessInfinity; ...
 'Degree loss stays separate from finite-root instability',@witnessDegreeLoss; ...
 'A stable report has no fabricated witness',@noWitness; ...
 'Diagnostic examples retain ordinary data and serialize without handles',@exampleData; ...
 'Diagnostic payload round-trips through the unchanged two-entry report ZIP',@()roundtrip(writer)};
end
function bernstein()
[~,E] = qstab.valueenclosure({[1 1],[3 1]},1i);r=E.regions{1};
assert(max(abs(r.controls-[2i;2+4i;8+6i]))<1e-12);
[~,S] = qstab.valueset({[1 1],[3 1]},1i,'Weights',[.5 .5]);assert(abs(S.points-(3+4i))<1e-12);
assert(~E.isCertificate&&~E.validated&&strcmp(E.parameterCoverage,'fullOriginalDomain'));
end
function falseAlarm()
V=alarmFamily();[~,E] = qstab.valueenclosure(V,0);r=E.regions{1};
assert(~E.zeroExcludedNumerically&&r.distanceToZero==0);
assert(min(real(r.controls))==-3&&max(real(r.controls))==5);
[~,S] = qstab.valueset(V,0,'Weights',[0 1;.5 .5;1 0]);assert(all(real(S.points)>=1));
end
function splitHull()
V=alarmFamily();cells={[1,.5;0,.5],[.5,0;.5,1]};
[~,E] = qstab.valueenclosure(V,0,'VertexCells',cells);
assert(numel(E.regions)==2&&E.zeroExcludedNumerically);
for j=1:2,assert(E.regions{j}.distanceToZero>=1-1e-12);assert(isequal(E.regions{j}.vertexWeights,cells{j}));end
end
function coverageFlag()
[~,E] = qstab.valueenclosure({[1 1],[2 1]},1i,'VertexCells',{[1;0]});
assert(strcmp(E.parameterCoverage,'userSubcellsCoverageNotChecked')&&~E.isCertificate);
end
function degenerateHull()
[~,E] = qstab.valueenclosure({3},1i);assert(isequal(E.regions{1}.polygon,9));
[~,E] = qstab.valueenclosure({1,2},0);assert(numel(E.regions{1}.polygon)==2);
[~,E] = qstab.valueenclosure({0},0);assert(~E.zeroExcludedNumerically);
end
function evenOdd()
[~,E] = qstab.valueenclosure({[1 2],[3 4]},1i,'Uncertainty','interval');r=E.regions{1};
assert(strcmp(r.kind,'EvenOddBoundingRectangle'));
assert(isequal(r.valueBounds,[-15 5;4 24]));
assert(E.zeroExcludedNumerically);
end
function correlated()
F=struct('center',poly([0 0;1 1;0 0;0 0]),'generators',{{[1 1]}},'bounds',[-2 2]);
[~,E] = qstab.valueenclosure(F,1i);assert(max(abs(E.regions{1}.valueBounds(:)-[0;2;0;10]))<1e-12);
assert(E.zeroExcludedNumerically);
end
function unitCircle()
z=exp(.7i);[~,E] = qstab.valueenclosure({[-.3 1],[-.1 1]},z,'Uncertainty','interval');
assert(strcmp(E.regions{1}.kind,'JointQuadraticRectangle'));
for t=linspace(0,1,21)
    q=(z-.3+.2*t)^2;insideRectangle(q,E.regions{1}.valueBounds);
end
end
function inputVariable()
[~,E] = qstab.valueenclosure({[1 -.5]},0);assert(E.regions{1}.controls==1);
end
function parameterBoxes()
[~,E] = qstab.valueenclosure({[1 1],[2 1]},1i,'Uncertainty','interval', ...
    'ParameterBoxes',{[1 1.5],[1.5 2]});
assert(numel(E.regions)==2&&isequal(E.regions{1}.parameterBounds,[1 1.5]));
assert(strcmp(E.parameterCoverage,'userSubcellsCoverageNotChecked'));
end
function sides()
p0=poly([2 1;1 0;0 0;0 0]);p1=poly([0;0;1;0]);
B=[0 0;1 1;0 0;0 0];F=struct('p0',p0,'p1',p1,'side','left','bounds',B);
F.p0=poly([2 1;0 0;0 0;1 0]);
[~,L] = qstab.valueenclosure(F,1i);F.side='right';[~,R] = qstab.valueenclosure(F,1i);
assert(norm(L.regions{1}.polygon-R.regions{1}.polygon)>1);
end
function invalidInputs()
expect(@()qstab.valueenclosure({[1 1]},Inf),'qstab:valueenclosure:Lambda');
expect(@()qstab.valueenclosure({[1 1]},[0 1]),'qstab:valueenclosure:Lambda');
expect(@()qstab.valueenclosure({[1 1],[2 1]},0,'VertexCells',{[-1;2]}),'qstab:valueenclosure:Weights');
expect(@()qstab.valueenclosure({[1 1],[2 1]},0,'VertexCells',{[1;1]}),'qstab:valueenclosure:Weights');
expect(@()qstab.valueenclosure({[1 1],[2 1]},0,'Uncertainty','interval', ...
    'ParameterBoxes',{[0 3]}),'qstab:valueenclosure:Bounds');
end
function overflow()
expect(@()qstab.valueenclosure({[1e300 1]},1i),'qstab:valueenclosure:Overflow');
end
function noSideEffects()
V=alarmFamily();old=V;state=rng;directory=pwd;
qstab.valueenclosure(V,1i);assert(isequal(V,old)&&isequal(state,rng)&&strcmp(directory,pwd));
end
function retained()
R=stablePoly();old=R;[~,D] = qstab.certificateData(R);
assert(D.reportedComplete&&D.storageComplete&&D.storedAcceptedCount==R.certifiedCells);
assert(strcmp(D.status,R.status)&&isequaln(R,old));
end
function missingStorage()
[~,R] = qstab.zeroexclude({[1 1],[2 1]},'MaxTime',Inf);[~,D] = qstab.certificateData(R);
assert(strcmp(D.status,'stable')&&D.reportedComplete&&~D.storageComplete&&isempty(D.rows));
end
function pending()
[~,R] = qstab.zeroexclude(repmat({[1 1]},1,5),'Reduction','four-generator', ...
    'MaxSubfamilies',1,'StoreCertificates',true,'MaxTime',Inf);
[~,D] = qstab.certificateData(R);assert(strcmp(D.status,'undetermined')&&~D.reportedComplete);
assert(any(strcmp({D.rows.state},'unresolved')));
assert(any(arrayfun(@(r)strcmp(r.state,'unresolved')&&isempty(r.frequencyBounds),D.rows)));
end
function lyapunov()
[~,R] = qstab.robuststable({[1 1],[1.1 1]},'StoreCertificates',true,'MaxTime',Inf);
[~,D] = qstab.certificateData(R);assert(D.storageComplete);
assert(all(strcmp({D.rows.metricUnits},'normalizedLyapunovInequality')));
assert(all(arrayfun(@(r)isempty(r.frequencyBounds),D.rows)));
end
function constant()
[~,R] = qstab.zeroexclude({1,2},'StoreCertificates',true,'MaxTime',Inf);[~,D] = qstab.certificateData(R);
assert(strcmp(D.status,'stable')&&isempty(D.rows)&&~D.storageComplete);
end
function certificatePoly()
V={[1 1],[2 1]};R=stablePoly();[~,E] = qstab.certificateEnclosure(V,R);
assert(~E.isCertificate&&strcmp(E.scope,'oneStoredParameterFrequencyCell'));
assert(E.geometryResidual<1e-7&&E.regions{1}.inflationRadius>=0);
end
function certificateFour()
V=repmat({[1 1]},1,5);[~,R] = qstab.zeroexclude(V,'Reduction','four-generator', ...
    'StoreCertificates',true,'MaxTime',Inf);
[~,D] = qstab.certificateData(R);assert(numel(D.rows(1).vertexWeights(:,1))==5);
[~,E] = qstab.certificateEnclosure(V,R);assert(E.geometryResidual<1e-7);
end
function wrongFamily()
R=stablePoly();expect(@()qstab.certificateEnclosure({[10 10],[20 10]},R),'qstab:certificateEnclosure:FamilyMismatch');
end
function certificateRectangles()
[~,R] = qstab.zeroexclude({[1 1],[1.1 1]},'Uncertainty','interval', ...
    'StoreCertificates',true,'MaxTime',Inf);
[~,E] = qstab.certificateEnclosure([],R);assert(strcmp(E.inputBinding,'geometryReadFromReport'));
F=struct('center',[1 1],'generators',{{.1}},'bounds',[-1 1]);
[~,R] = qstab.zeroexclude(F,'StoreCertificates',true,'MaxTime',Inf);
[~,E] = qstab.certificateEnclosure([],R);assert(E.regions{1}.inflationRadius>=0);
end
function certificateDelay()
V={[1 -.2],[1 -.3]};[~,R] = qstab.zeroexclude(V,'Variable','zinv', ...
    'StoreCertificates',true,'MaxTime',Inf);[~,E] = qstab.certificateEnclosure(V,R);
assert(strcmp(E.originalVariable,'zinv')&&strcmp(E.variable,'z'));
assert(abs(abs(E.lambda)-1)<1e-12&&E.geometryResidual<1e-7);
end
function outerDisk()
p=qstab.internal.diagnosticGeometry('inflate',0,1,32);
% Polygon sides are tangent to the unit circle; vertices lie outside it.
assert(all(abs(p)>1));
for j=1:numel(p)
    a=p(j);b=p(mod(j,numel(p))+1);mid=(a+b)/2;assert(abs(abs(mid)-1)<1e-12);
end
end
function noGeometry()
[~,R] = qstab.robuststable({[1 1]},'StoreCertificates',true,'MaxTime',Inf);
expect(@()qstab.certificateEnclosure([],R),'qstab:certificateEnclosure:NoGeometry');
R=stablePoly();R.certificates{1}.cell.frequencyBounds=[0 1 2];
expect(@()qstab.certificateEnclosure({[1 1],[2 1]},R),'qstab:certificateEnclosure:Frequency');
end
function witnessPoly()
V={[-1 1],[-2 1]};[~,R] = qstab.zeroexclude(V,'MaxTime',Inf);[~,W] = qstab.witnessData(V,R);
assert(W.isAdmissible&&W.isConfirmedCounterexample&&W.spectrumValid);
assert(abs(W.inputRepresentatives-1.5)<1e-10&&max(abs(W.weights-.5))<1e-14);
end
function witnessTampered()
V={[-1 1],[-2 1]};[~,R] = qstab.zeroexclude(V,'MaxTime',Inf);
R.witness.components(1,1)=-5;[~,W] = qstab.witnessData(V,R);assert(~W.isAdmissible&&isempty(W.nominal));
[~,R] = qstab.zeroexclude(V,'MaxTime',Inf);R.witness.originalWeights=[-1;2];
[~,W] = qstab.witnessData(V,R);assert(~W.isAdmissible);
end
function witnessInterval()
V={[-2 1],[-1 1]};[~,R] = qstab.zeroexclude(V,'Uncertainty','interval','MaxTime',Inf);
[~,W] = qstab.witnessData(V,R);assert(W.isAdmissible&&W.isConfirmedCounterexample);
R.witness.components(2,1)=1;[~,W] = qstab.witnessData(V,R);assert(~W.isAdmissible);
end
function witnessAffine()
F=struct('center',[-1 1],'generators',{{[.1 0]}},'bounds',[-1 1]);
[~,R] = qstab.zeroexclude(F,'MaxTime',Inf);[~,W] = qstab.witnessData(F,R);
assert(W.isAdmissible&&W.isConfirmedCounterexample&&W.parameters==0);
end
function witnessDelay()
V={[1 -1.2],[1 -1.4]};[~,R] = qstab.zeroexclude(V,'Variable','zinv','MaxTime',Inf);
[~,W] = qstab.witnessData(V,R);assert(W.isConfirmedCounterexample);
assert(abs(W.testRepresentatives-1.3)<1e-10&&abs(W.inputRepresentatives-1/1.3)<1e-10);
assert(strcmp(W.variable,'zinv')&&strcmp(W.testVariable,'z'));
end
function witnessInfinity()
V={[1 -1.2 0],[1 -1.4 0]};[~,R] = qstab.zeroexclude(V,'Variable','zinv','MaxTime',Inf);
[~,W] = qstab.witnessData(V,R);assert(W.isConfirmedCounterexample&&W.grade==2);
assert(sum(isinf(abs(W.inputRepresentatives)))==1&&sum(W.testRepresentatives==0)==1);
end
function witnessDegreeLoss()
V={[1 1],[-1 -1]};[~,R] = qstab.zeroexclude(V,'MaxTime',Inf);[~,W] = qstab.witnessData(V,R);
assert(W.isAdmissible&&W.isDegreeLoss&&strcmp(W.checkStatus,'confirmedDegreeLoss'));
assert(~W.isConfirmedCounterexample&&isempty(W.inputRepresentatives));
end
function noWitness()
R=stablePoly();[~,W] = qstab.witnessData({[1 1],[2 1]},R);
assert(~W.available&&~W.isAdmissible&&strcmp(W.checkStatus,'noWitness'));
end
function exampleData()
D=run_qstab_diagnostics('Plots',false,'Verbose',false);assert(D.ok);
assert(~containsHandle(D)&&strcmp(D.schema,'qstab.diagnostic-examples.1'));
end
function roundtrip(writer)
base=tempname;mkdir(base);cleanup=onCleanup(@()rmdir(base,'s')); %#ok<NASGU>
item=struct('group','fixture','name','Artificial diagnostic persistence fixture', ...
    'status','PASS','seconds',0,'message','','identifier','');
r=struct('toolbox',struct('version','fixture'), ...
    'environment',struct('version','artificial','computer','fixture','quaternionPath','','qstabPath',''), ...
    'timestamp','20000101T000000','runId','diagnostic-persistence-fixture', ...
    'nativeMATLABRun',false,'nativeQuaternionOK',false,'ok',true,'status','OVERALL_OK', ...
    'passed',1,'failed',0,'skipped',0,'elapsedSeconds',0,'cases',item);
r.diagnosticExamples=run_qstab_diagnostics('Plots',false,'Verbose',false);
r=writer(r,base,'qstab_test_results');
assert(strcmp(r.archiveFile,fullfile(base,'qstab_test_results.zip')));
folder=fullfile(base,'unpacked');mkdir(folder);unzip(r.archiveFile,folder);
d=dir(folder);d=d(~ismember({d.name},{'.','..'}));
assert(numel(d)==2&&isequal(sort({d.name}),sort({'report.mat','report.txt'})));
s=load(fullfile(folder,'report.mat'),'report');assert(isequaln(s.report.diagnosticExamples,r.diagnosticExamples));
end
function R=stablePoly()
[~,R] = qstab.zeroexclude({[1 1],[2 1]},'StoreCertificates',true,'MaxTime',Inf);assert(strcmp(R.status,'stable'));
end
function V=alarmFamily()
V={poly([1 1;0 0;-2 0;0 0]),poly([1 1;0 0;2 0;0 0])};
end
function P=poly(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
function insideRectangle(z,B)
t=1e-11*max(1,abs(z));assert(real(z)>=B(1,1)-t&&real(z)<=B(1,2)+t&&imag(z)>=B(2,1)-t&&imag(z)<=B(2,2)+t);
end
function expect(f,id)
try,f();catch ME,assert(strcmp(ME.identifier,id),sprintf('Expected %s, got %s',id,ME.identifier));return;end
error('qstab:test:MissingError','Expected error %s.',id);
end
function tf=containsHandle(x)
tf=isa(x,'handle');if tf,return;end
if isstruct(x)
    names=fieldnames(x);for a=1:numel(x),for b=1:numel(names),if containsHandle(x(a).(names{b})),tf=true;return;end;end;end
elseif iscell(x)
    for a=1:numel(x),if containsHandle(x{a}),tf=true;return;end;end
end
end
