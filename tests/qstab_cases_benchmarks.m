function cases = qstab_cases_benchmarks(writer)
%QSTAB_CASES_BENCHMARKS Benchmark policy, not an additional stability theorem.
% Artificial timing/status records below test reporting rules only. They
% are never included in the production benchmark as measured solver data.
cases={ ...
 'Benchmark fixtures have unique IDs, explicit references and no random inputs',@fixtureContract; ...
 'Every method pairing is applicable or explicitly rejected without fallback',@applicability; ...
 'All applicable rows receive identical common tolerances and budgets',@matchedOptions; ...
 'CT cube and quadratic restrictions are preserved by the benchmark plan',@restrictions; ...
 'Delay fixtures are exact coefficient reversals at unchanged grade',@delayInputs; ...
 'One-quaternion inputs preserve their multiplication side',@sideInputs; ...
 'A stable benchmark result requires complete reported coverage',@stableCoverage; ...
 'Unresolved results are retained and never converted to instability',@unknown; ...
 'Resolved statuses contradicting an analytic reference are rejected',@contradiction; ...
 'A polytope counterexample reconstructs all original convex weights',@polyWitness; ...
 'Out-of-family or mismatched counterexamples are rejected',@badWitness; ...
 'Delay witnesses are rechecked in the original input variable',@delayWitness; ...
 'An interval witness is checked against every component bound',@intervalWitness; ...
 'Degree loss remains separate from instability and timing rankings',@degreeWitness; ...
 'Solved timing ratios exclude unknown and notApplicable rows',@resolvedRatios; ...
 'Short timing batches do not produce a speedup claim',@shortTiming; ...
 'A warm-up error cannot disappear from the benchmark summary',@warmupFailure; ...
 'Mixed budget verdicts remain visible and unranked',@mixedVerdicts; ...
 'Conflicting resolved verdicts flag a benchmark issue',@conflictingVerdicts; ...
 'A small real benchmark runs both methods and retains complete outputs',@smoke; ...
 'Method order rotates without reseeding the random number generator',@orderAndRng; ...
 'A subfamily budget is recorded as unknown rather than a faster solution',@bounded; ...
 'An entirely unsupported selection performs no solver calls',@unsupported; ...
 'Malformed benchmark options are rejected',@invalid; ...
 'A return-only benchmark does not create archives or change the directory',@noIO; ...
 'Benchmark payload round-trips through the existing two-entry report archive',@()persistence(writer)};
end
function o=options()
o=struct('Methods',{{}},'Tolerance',1e-10,'MaxCells',512, ...
 'MaxDepth',32,'MaxTime',1,'MaxSubfamilies',32,'Repetitions',1,'MinSampleTime',.1);
end
function c=fixture(id)
a=qstab_benchmark_cases();k=find(strcmp({a.id},id));assert(numel(k)==1);c=a(k);
end
function fixtureContract()
a=qstab_benchmark_cases();assert(numel(a)==16&&numel(unique({a.id}))==16);
for j=1:numel(a)
 F=qstab.familycomponents(a(j).data,'Uncertainty',a(j).uncertainty);
 assert(F.degree>=1&&~isempty(a(j).referenceReason)&&~isempty(a(j).source));
end
assert(isequaln(a,qstab_benchmark_cases()));
end
function applicability()
a=qstab_benchmark_cases();p=qstab_benchmark_plan(a,options());
assert(numel(p)==96&&sum([p.applicable])==49);
assert(all(cellfun(@(s)~isempty(s),{p.reason})));
end
function matchedOptions()
p=qstab_benchmark_plan(qstab_benchmark_cases(),options());o=options();
fields={'Tolerance','MaxCells','MaxDepth','MaxTime','MaxSubfamilies'};
for j=1:numel(p)
 for k=1:numel(fields),assert(isequal(p(j).solverOptions.(fields{k}),o.(fields{k})));end
 assert(~p(j).solverOptions.StoreCertificates);
end
end
function restrictions()
p=qstab_benchmark_plan(qstab_benchmark_cases(),options());
for j=1:numel(p)
 if strcmp(p(j).method,'quadratic-faces')&&p(j).applicable
  assert(any(strcmp(p(j).caseId,{'quadratic_interval','quadratic_interior'})));
 end
 if strncmp(p(j).method,'kharitonov',9)&&p(j).applicable
  assert(strcmp(p(j).solverOptions.Variable,'s')&&strcmp(p(j).solverOptions.Uncertainty,'interval'));
 end
end
end
function delayInputs()
a=fixture('schur_polytope5');b=fixture('delay_polytope5');
A=qstab.familycomponents(a.data);B=qstab.familycomponents(b.data);
assert(isequal(A.vertices(:,end:-1:1,:),B.vertices));
a=fixture('schur_interval');b=fixture('delay_interval');
A=qstab.familycomponents(a.data,'Uncertainty','interval');B=qstab.familycomponents(b.data,'Uncertainty','interval');
assert(isequal(A.intervalLower(:,end:-1:1),B.intervalLower));
assert(isequal(A.intervalUpper(:,end:-1:1),B.intervalUpper));
end
function sideInputs()
a=fixture('oneparameter_left');b=fixture('oneparameter_right');
assert(strcmp(a.data.side,'left')&&strcmp(b.data.side,'right'));
A=qstab.familycomponents(a.data);B=qstab.familycomponents(b.data);
assert(~isequal(A.generators,B.generators)&&isequal(A.bounds,B.bounds));
end
function R=result(status)
R=struct('status',status,'isStable',[],'certificateComplete',false, ...
 'unresolvedCells',0,'method','syntheticPolicyFixture','witness',[]);
if strcmp(status,'stable'),R.isStable=true;R.certificateComplete=true;
elseif strcmp(status,'notStable'),R.isStable=false;end
end
function stableCoverage()
c=fixture('paired_interval');r=result('stable');a=qstab_benchmark_assess(c,r,1e-10);
assert(a.ok&&a.comparisonEligible);
r.certificateComplete=false;a=qstab_benchmark_assess(c,r,1e-10);assert(~a.ok);
r.certificateComplete=true;r.unresolvedCells=1;a=qstab_benchmark_assess(c,r,1e-10);assert(~a.ok);
end
function unknown()
c=fixture('paired_interval');r=result('undetermined');a=qstab_benchmark_assess(c,r,1e-10);
assert(a.ok&&~a.comparisonEligible&&~a.resolved);
r.isStable=false;a=qstab_benchmark_assess(c,r,1e-10);assert(~a.ok);
end
function contradiction()
r=result('stable');a=qstab_benchmark_assess(fixture('quadratic_interior'),r,1e-10);assert(~a.ok);
a=qstab_benchmark_assess(fixture('boundary_touch'),r,1e-10);assert(~a.ok);
a=qstab_benchmark_assess(fixture('degree_loss'),r,1e-10);assert(~a.ok);
end
function [c,r]=negativePoly()
c=fixture('tetrahedron5');F=qstab.familycomponents(c.data);w=[.25;.25;.25;.25;0];
C=zeros(4,3);for j=1:5,C=C+w(j)*F.vertices(:,:,j);end
r=result('notStable');r.witness=struct('originalComponents',C,'originalWeights',w);
end
function polyWitness()
[c,r]=negativePoly();a=qstab_benchmark_assess(c,r,1e-10);
assert(a.ok&&a.witnessMember&&a.comparisonEligible&&strcmp(a.witnessNominalStatus,'notStable'));
end
function badWitness()
[c,r]=negativePoly();r.witness.originalWeights(1)=-.25;
a=qstab_benchmark_assess(c,r,1e-10);assert(~a.ok);
[c,r]=negativePoly();r.witness.originalComponents(1,1)=100;
a=qstab_benchmark_assess(c,r,1e-10);assert(~a.ok);
end
function delayWitness()
c=fixture('delay_polytope5');c.data={[1,2]};c.referenceStatus='notStable';
r=result('notStable');r.witness=struct('components',[2,1;0,0;0,0;0,0],'originalWeights',1);
a=qstab_benchmark_assess(c,r,1e-10);assert(a.ok&&a.witnessMember);
end
function intervalWitness()
c=fixture('quadratic_interior');F=qstab.familycomponents(c.data,'Uncertainty','interval');
r=result('notStable');r.witness=struct('originalComponents',.5*F.intervalLower+.5*F.intervalUpper);
a=qstab_benchmark_assess(c,r,1e-10);assert(a.ok&&a.witnessMember);
r.witness.originalComponents(3,2)=3;
a=qstab_benchmark_assess(c,r,1e-10);assert(~a.ok);
end
function degreeWitness()
c=fixture('degree_loss');r=result('assumptionFailed');
r.witness=struct('originalComponents',zeros(4,2),'originalWeights',[.5;.5]);
a=qstab_benchmark_assess(c,r,1e-10);assert(a.ok&&~a.comparisonEligible);
end
function x=observation(status,eligible,ok,t)
x=struct('status',status,'comparisonEligible',eligible,'ok',ok, ...
 'seconds',t,'processedCells',10,'processedSubfamilies',1, ...
 'result',struct('purpose','synthetic benchmark-policy test, not measured evidence'));
end
function r=row(obs,time)
s=struct('calls',numel(obs),'elapsedSeconds',time, ...
 'secondsPerCall',time/numel(obs),'observations',{obs});
r=struct('warmups',{{}},'samples',{{s}});
end
function [p,rows,o]=summaryFixture()
o=options();o.Methods={'direct','lyapunov'};p=qstab_benchmark_plan(fixture('paired_interval'),o);
rows={row({observation('stable',true,true,.2)},.2); ...
      row({observation('stable',true,true,.4)},.4)};
end
function resolvedRatios()
[p,r,o]=summaryFixture();[s,q]=qstab_benchmark_summarize(p,r,o);
assert(all([s.comparisonEligible])&&q.eligible&&abs(q.timeAoverB-.5)<eps);
r{2}=row({observation('undetermined',false,true,.01)},.01);
[s,q]=qstab_benchmark_summarize(p,r,o);assert(~q.eligible&&isnan(q.timeAoverB)&&s(2).ok);
end
function shortTiming()
[p,r,o]=summaryFixture();r{1}=row({observation('stable',true,true,.01)},.01);
[s,q]=qstab_benchmark_summarize(p,r,o);assert(s(1).decided&&~s(1).timingAdequate&&~q.eligible);
end
function warmupFailure()
[p,r,o]=summaryFixture();r{1}.warmups={observation('error',false,false,.01)};
[s,q]=qstab_benchmark_summarize(p,r,o);assert(~s(1).ok&&~q.eligible);
end
function mixedVerdicts()
[p,r,o]=summaryFixture();r{1}=row({observation('stable',true,true,.1); ...
 observation('undetermined',false,true,.1)},.2);
[s,q]=qstab_benchmark_summarize(p,r,o);assert(s(1).ok&&strcmp(s(1).status,'mixed')&&~q.eligible);
end
function conflictingVerdicts()
[p,r,o]=summaryFixture();r{1}=row({observation('stable',true,true,.1); ...
 observation('notStable',true,true,.1)},.2);
[s,q]=qstab_benchmark_summarize(p,r,o);assert(~s(1).ok&&~q.eligible);
end
function B=small(varargin)
B=run_qstab_benchmarks('Cases',{'paired_interval'},'Methods',{'direct','lyapunov'}, ...
 'WarmupRuns',0,'Repetitions',1,'MinSampleTime',0,'MaxBatchCalls',1, ...
 'MaxTime',Inf,'Verbose',false,varargin{:});
end
function smoke()
B=small();assert(B.ok&&B.allDecided&&B.measuredSolverCalls==2&&numel(B.rows)==2);
assert(all(strcmp({B.summary.status},'stable'))&&all([B.summary.comparisonEligible]));
for k=1:2,assert(isstruct(B.rows{k}.samples{1}.observations{1}.result));end
end
function orderAndRng()
state=rng;B=small('Repetitions',3);after=rng;
assert(isequal(state,after)&&B.measuredSolverCalls==6);
assert(isequal(B.schedule(:,2:3),[1,1;2,1;2,2;1,2;1,3;2,3]));
end
function bounded()
B=small('Cases',{'polytope_ct6'},'Methods',{'four-generator'},'MaxSubfamilies',1);
assert(B.ok&&~B.allDecided&&strcmp(B.summary.status,'undetermined')&&~B.summary.comparisonEligible);
end
function unsupported()
B=small('Cases',{'schur_interval'},'Methods',{'quadratic-faces'});
assert(B.ok&&~B.allDecided&&B.measuredSolverCalls==0&&strcmp(B.summary.status,'notApplicable'));
end
function invalid()
bad={{'Cases',{'absent'}},{'Methods',{'bad'}},{'Cases',{'paired_interval','paired_interval'}}, ...
 {'Repetitions',0},{'MinSampleTime',-1},{'MaxBatchCalls',0},{'MaxCells',0}, ...
 {'Tolerance',NaN},{'WarmupRuns',-.1},{'Methods','direct'}};
for k=1:numel(bad)
 args=bad{k};caught=false;try,small(args{:});catch,caught=true;end
 assert(caught,'Malformed benchmark selector/budget must fail.');
end
end
function noIO()
old=pwd;root=tempname;mkdir(root);cleanup=onCleanup(@()tidy(old,root)); %#ok<NASGU>
cd(root);here=pwd;before=dir(root);before={before.name};B=small();after=dir(root);after={after.name};
assert(B.ok&&strcmp(pwd,here)&&isequal(before,after));
end
function persistence(writer)
base=tempname;mkdir(base);cleanup=onCleanup(@()rmdir(base,'s')); %#ok<NASGU>
item=struct('group','fixture','name','Artificial persistence fixture', ...
 'status','PASS','seconds',0,'message','','identifier','');
r=struct('toolbox',struct('version','fixture'), ...
 'environment',struct('version','artificial','computer','fixture','quaternionPath','','qstabPath',''), ...
 'timestamp','20000101T000000','runId','benchmark-persistence-fixture', ...
 'nativeMATLABRun',false,'nativeQuaternionOK',false,'ok',true,'status','OVERALL_OK', ...
 'passed',1,'failed',0,'skipped',0,'elapsedSeconds',0,'cases',item);
r.benchmarks=small();r=writer(r,base,'qstab_test_results');
assert(strcmp(r.archiveFile,fullfile(base,'qstab_test_results.zip')));
folder=fullfile(base,'unpacked');mkdir(folder);unzip(r.archiveFile,folder);
d=dir(folder);d=d(~ismember({d.name},{'.','..'}));assert(numel(d)==2);
s=load(fullfile(folder,'report.mat'),'report');assert(isequaln(s.report.benchmarks,r.benchmarks));
end
function tidy(old,root)
cd(old);rmdir(root,'s');
end
