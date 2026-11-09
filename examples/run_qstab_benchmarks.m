function B = run_qstab_benchmarks(varargin)
%RUN_QSTAB_BENCHMARKS Matched, bounded scalar robustness comparison.
% B=run_qstab_benchmarks returns ordinary structs only; it does NOT write
% files, rename archives, change paths, seed RNG, or time plot/export work.
% To save it in the usual two-entry archive use:
%   report=test_qstab('Benchmarks',true);
%   % report.benchmarks is inside qstab_public/qstab_test_results.zip.
%
% Methods: lyapunov, direct, four-generator, kharitonov-cube,
% kharitonov-faces, quadratic-faces. Unsupported pairs are recorded as
% notApplicable, never silently routed to a different method.
% Same family and Tolerance/MaxCells/MaxDepth/MaxTime/MaxSubfamilies apply
% to all compared calls. Cells from different engines are NOT equal work.
% Full-family reductions retain their native GLOBAL subproblem budgets.
%
% Defaults: 3 repetitions, 1 warm-up, MinSampleTime=.1, MaxBatchCalls=32,
% MaxCells=512, MaxDepth=32, MaxTime=1, MaxSubfamilies=32, Tolerance=1e-10.
% Explicit tic handles are used because the solvers contain their own
% tic/toc budget timers. Every timed output is consumed and retained;
% witness checking is outside measured time. Method order is rotated.
% A short batch is labelled timingAdequate=false, not silently accepted.
% No universal faster-method claim or cross-machine comparison is made.
%
% Cases/Methods select cell vectors of exact IDs. List cases using
% qstab_benchmark_cases. Warmups and raw timed results are kept in B.rows.
% Unknown/budget-limited is valid data; runtime errors or contradictions
% set B.ok=false. Completion does NOT imply every family was decided.
started=tic;ip=inputParser;ip.FunctionName='run_qstab_benchmarks';
addParameter(ip,'Cases',{},@textList);addParameter(ip,'Methods',{},@textList);
addParameter(ip,'Repetitions',3,@(x)integer(x,1));
addParameter(ip,'WarmupRuns',1,@(x)integer(x,0));
addParameter(ip,'MinSampleTime',.1,@(x)finiteScalar(x,0,Inf));
addParameter(ip,'MaxBatchCalls',32,@(x)integer(x,1));
addParameter(ip,'Tolerance',1e-10,@(x)finiteScalar(x,realmin,.01)&&x<.01);
addParameter(ip,'MaxCells',512,@(x)integer(x,1));
addParameter(ip,'MaxDepth',32,@(x)integer(x,0));
addParameter(ip,'MaxTime',1,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&x>=0&&~isnan(x));
addParameter(ip,'MaxSubfamilies',32,@(x)integer(x,1));
addParameter(ip,'Verbose',true,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});o=ip.Results;
o.Cases=normalize(o.Cases);o.Methods=normalize(o.Methods);
allCases=qstab_benchmark_cases();
if isempty(o.Cases),cases=allCases;
else
    ids={allCases.id};[found,indices]=ismember(o.Cases,ids);
    if ~all(found),error('qstab:benchmark:Case','Unknown benchmark case.');end
    cases=allCases(indices);
end
plan=qstab_benchmark_plan(cases,o);
B=struct('schema','qstab.benchmark.1','toolbox',qstab.about(), ...
    'environment',struct('version',version,'computer',computer, ...
        'isOctave',exist('OCTAVE_VERSION','builtin')~=0), ...
    'options',o,'cases',cases,'plan',plan,'rows',{{}}, ...
    'summary',[],'comparisons',[],'ok',false,'allDecided',false, ...
    'status','RUNNING','elapsedSeconds',0,'measuredSolverCalls',0, ...
    'warmupSolverCalls',0,'schedule',zeros(0,3),'sourceRunId','', ...
    'assurance','ordinary double; numerical reference/witness checks, not interval proofs', ...
    'timingPolicy','median of explicit-handle wall-time batches after warm-up; order rotated', ...
    'archivePolicy','return-only; test_qstab attaches B to report.mat in the unchanged ZIP');
B.rows=cell(numel(plan),1);
for j=1:numel(plan)
    B.rows{j}=struct('warmups',{{}},'samples',{{}});
end
serial=0;
for c=1:numel(cases)
    indices=find([plan.caseIndex]==c & [plan.applicable]);m=numel(indices);
    if m==0,continue;end
    % Warm-up uses the same one-output call and the same numerical limits.
    for w=1:o.WarmupRuns
        for k=1:m
            j=indices(k);obs=call(cases(c),plan(j),o.Tolerance);
            B.rows{j}.warmups{end+1,1}=obs;B.warmupSolverCalls=B.warmupSolverCalls+1;
        end
    end
    for rep=1:o.Repetitions
        order=indices(mod((0:m-1)+(rep-1),m)+1);
        for j=reshape(order,1,[])
            elapsed=0;obs={};calls=0;
            while calls<o.MaxBatchCalls && (calls==0||elapsed<o.MinSampleTime)
                entry=call(cases(c),plan(j),o.Tolerance);
                calls=calls+1;elapsed=elapsed+entry.seconds;obs{end+1,1}=entry; %#ok<AGROW>
                if ~entry.ok,break;end % keep the failed observation, no retry-until-pass
            end
            sample=struct('repetition',rep,'calls',calls,'elapsedSeconds',elapsed, ...
                'secondsPerCall',elapsed/calls,'observations',{obs});
            B.rows{j}.samples{end+1,1}=sample;
            B.measuredSolverCalls=B.measuredSolverCalls+calls;
            serial=serial+1;B.schedule(end+1,:)=[serial,j,rep]; %#ok<AGROW>
        end
    end
    if o.Verbose,fprintf('Benchmark fixture completed: %s\n',cases(c).id);end
end
[B.summary,B.comparisons]=qstab_benchmark_summarize(plan,B.rows,o);
applicable=[B.summary.applicable];B.ok=all([B.summary(applicable).ok]);
B.allDecided=any(applicable)&&all([B.summary(applicable).decided]);
if ~any(applicable),B.status='NO_APPLICABLE_METHODS';
elseif B.ok,B.status='COMPLETED';else,B.status='COMPLETED_WITH_ISSUES';end
B.elapsedSeconds=toc(started);
if o.Verbose
    fprintf('\nqstab matched benchmark: %s\n',B.status);
    fprintf('%-24s %-18s %-16s %10s %8s %8s\n','family','method','status','median s','cells','calls');
    for j=1:numel(B.summary)
        s=B.summary(j);if ~s.applicable,continue;end
        fprintf('%-24s %-18s %-16s %10.5g %8g %8g\n', ...
            s.caseId,s.method,s.status,s.medianSeconds,s.medianCells,s.measuredCalls);
    end
    fprintf('Stored %d applicable and %d notApplicable rows; %d warm-up / %d measured calls.\n', ...
        sum(applicable),sum(~applicable),B.warmupSolverCalls,B.measuredSolverCalls);
    fprintf('Unresolved is not failure and is never ranked as a faster solution.\n');
    fprintf('Raw inputs, options, outputs, timing batches and admissible witnesses are retained.\n');
end
end
function obs=call(fixture,p,tolerance)
% No input conversion, witness verification, reporting or I/O is timed.
names=fieldnames(p.solverOptions);args=cell(1,2*numel(names));
for k=1:numel(names),args{2*k-1}=names{k};args{2*k}=p.solverOptions.(names{k});end
obs=struct('status','error','ok',false,'comparisonEligible',false, ...
    'seconds',NaN,'processedCells',NaN,'processedSubfamilies',NaN, ...
    'result',[],'assessment',[],'errorIdentifier','','errorMessage','');
timer=tic;
try
    [~,R] = qstab.robuststable(fixture.data,args{:});obs.seconds=toc(timer);
catch ME
    obs.seconds=toc(timer);obs.errorIdentifier=ME.identifier;obs.errorMessage=ME.message;return
end
obs.result=R;obs.status=R.status;
if isfield(R,'processedCells'),obs.processedCells=R.processedCells;end
if isfield(R,'reduction')&&isfield(R.reduction,'processedSubfamilies')
    obs.processedSubfamilies=R.reduction.processedSubfamilies;
end
try
    A=qstab_benchmark_assess(fixture,R,tolerance);
    obs.assessment=A;obs.ok=A.ok;obs.comparisonEligible=A.comparisonEligible;
catch ME
    obs.errorIdentifier=ME.identifier;obs.errorMessage=['Assessment: ',ME.message];
end
end
function tf=textList(x)
tf=iscell(x)&&(isempty(x)||isvector(x));
if tf
    for k=1:numel(x)
        if ~((ischar(x{k})&&isrow(x{k})&&~isempty(strtrim(x{k})))|| ...
                (isstring(x{k})&&isscalar(x{k})&&strlength(x{k})>0))
            tf=false;return
        end
    end
end
end
function out=normalize(x)
out=cellfun(@(s)lower(strtrim(char(s))),reshape(x,1,[]),'UniformOutput',false);
if numel(unique(out))~=numel(out),error('qstab:benchmark:Duplicate','Duplicate selectors are not permitted.');end
end
function tf=integer(x,minimum)
tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=minimum&&x==floor(x);
end
function tf=finiteScalar(x,lo,hi)
tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=lo&&x<=hi;
end
