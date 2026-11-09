function T=qstab_eigen_timing(first,second,A,reference,varargin)
%QSTAB_EIGEN_TIMING Matched ONE-output batches, with retained call results.
% Warm-up/calibration and measured calls are separately recorded. Every
% output is checked outside the measured batch. Ratios are suppressed when
% either method is incorrect, a batch is short, or any invocation errors.
if ~isa(first,'function_handle')||~isa(second,'function_handle')
    error('qstab:eigencomparison:Handle','Two explicit function handles are required.');
end
ip=inputParser;
addParameter(ip,'Tolerance',1e-7,@positive);
addParameter(ip,'Repetitions',3,@posint);
addParameter(ip,'WarmupRuns',1,@nonnegint);
addParameter(ip,'MinSampleTime',0.05,@nonnegative);
addParameter(ip,'MaxBatchCalls',256,@posint);
parse(ip,varargin{:});o=ip.Results;
T=struct('status','RUNNING','ok',true,'eligible',false, ...
    'contract','exactly one output from each function; same input object and common batch count', ...
    'options',o,'warmup',{{}},'calibration',{{}},'samples',{{}}, ...
    'firstMedianSeconds',NaN,'secondMedianSeconds',NaN,'secondOverFirst',NaN, ...
    'firstCalls',0,'secondCalls',0,'message','');
for k=1:o.WarmupRuns
    B=pair(first,second,A,1,reference,o.Tolerance,mod(k,2)==0);
    T.warmup{end+1}=B;T=counts(T,B);
    if ~B.ok,T.status='CALL_OR_OUTPUT_ISSUE';T.ok=false;T.message='A warm-up failed; no timing ratio.';return;end
end
N=1;stage=0;
while true
    stage=stage+1;
    B=pair(first,second,A,N,reference,o.Tolerance,mod(stage,2)==0);
    T.calibration{end+1}=B;T=counts(T,B);
    if ~B.ok,T.status='CALL_OR_OUTPUT_ISSUE';T.ok=false;T.message='Calibration failed; no timing ratio.';return;end
    if min(B.first.seconds,B.second.seconds)>=o.MinSampleTime||N>=o.MaxBatchCalls,break;end
    N=min(o.MaxBatchCalls,2*N);
end
firstTimes=zeros(1,o.Repetitions);secondTimes=firstTimes;qualified=true;
for k=1:o.Repetitions
    B=pair(first,second,A,N,reference,o.Tolerance,mod(k,2)==0);
    T.samples{end+1}=B;T=counts(T,B);
    if ~B.ok,T.status='CALL_OR_OUTPUT_ISSUE';T.ok=false;T.message='A measured call failed; no timing ratio.';return;end
    firstTimes(k)=B.first.seconds/N;secondTimes(k)=B.second.seconds/N;
    qualified=qualified&&B.first.seconds>=o.MinSampleTime&&B.second.seconds>=o.MinSampleTime;
end
T.firstMedianSeconds=median(firstTimes);T.secondMedianSeconds=median(secondTimes);
T.eligible=qualified&&T.firstMedianSeconds>0&&T.secondMedianSeconds>0;
if T.eligible
    T.status='MEASURED';T.secondOverFirst=T.secondMedianSeconds/T.firstMedianSeconds;
    T.message='Both complete multisets agree with the reference and every measured batch met the time floor.';
else
    T.status='SHORT_BATCH';T.message='Raw timings retained; no ratio because one or more batches were too short.';
end
end
function B=pair(first,second,A,N,reference,tol,reverse)
if reverse
    S=batch(second,A,N,reference,tol);F=batch(first,A,N,reference,tol);order={'second','first'};
else
    F=batch(first,A,N,reference,tol);S=batch(second,A,N,reference,tol);order={'first','second'};
end
B=struct('callsPerMethod',N,'order',{order},'first',F,'second',S,'ok',F.ok&&S.ok);
end
function B=batch(f,A,N,reference,tol)
raw=cell(1,N);done=0;message='';identifier='';
t=tic;
try
    for k=1:N
        raw{k}=f(A); % Assigning is essential: nargout must be ONE, never zero.
        done=k;
    end
catch ME
    message=ME.message;identifier=ME.identifier;
end
elapsed=toc(t);
B=struct('seconds',elapsed,'requestedCalls',N,'completedCalls',done, ...
    'ok',isempty(message)&&done==N,'classes',{cell(1,done)},'errors',Inf(1,done), ...
    'message',message,'identifier',identifier);
for k=1:done
    try
        [a,b,c,d]=qstab.components(raw{k});
        if ~isvector(a)||numel(a)~=numel(reference)
            error('qstab:eigencomparison:TimingShape','A timed call returned the wrong eigenvalue shape.');
        end
        z=complex(a(:),hypot(hypot(b(:),c(:)),d(:)));
        M=qstab_eigen_match(z,reference,tol);B.classes{k}=z;B.errors(k)=M.error;
        B.ok=B.ok&&M.ok;
        if ~M.ok&&isempty(B.message),B.message=M.message;end
    catch ME
        B.ok=false;
        if isempty(B.message),B.message=ME.message;B.identifier=ME.identifier;end
    end
end
end
function T=counts(T,B)
T.firstCalls=T.firstCalls+B.first.completedCalls;
T.secondCalls=T.secondCalls+B.second.completedCalls;
end
function tf=positive(x),tf=nonnegative(x)&&x>0;end
function tf=posint(x),tf=positive(x)&&x==floor(x);end
function tf=nonnegint(x),tf=nonnegative(x)&&x==floor(x);end
function tf=nonnegative(x),tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=0;end
