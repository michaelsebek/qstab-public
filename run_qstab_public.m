function report=run_qstab_public(varargin)
%RUN_QSTAB_PUBLIC One standalone acceptance run, without an external qeig.
%   REPORT=run_qstab_public temporarily restricts the MATLAB search path to
%   this package and existing directories under MATLABROOT, runs installation
%   checks, the full regression suite and analytic eigenpair references, then
%   restores the caller's path and current folder (also on an exception).
%
%   Run from a freshly restarted PUBLIC MATLAB session. This routine does not
%   clear loaded classes/functions, reset persistent/global state or remove
%   user data. Path restriction is not a security/process sandbox.
%
%   Options (exact names, logical scalars):
%     Examples        true: all existing optional correctness/example groups
%                     except the long timing benchmark and figure export
%     Figures         false: additionally regenerate/export the eight figures
%                     and embed the PNG/PDF bytes in report.publicAcceptance
%     DependencyScan  true: MATLAB static dependencies of public src entries
%     Verbose         true: print ordinary regression/example progress
%
%   No LegacyFunction/EigenOptions override is accepted. qeig is never requested.
%   Timing is false; CheckPrimaryVectors and UseQuaternion are true. Native
%   quaternion and graphics regressions are required for acceptance.
%
%   The ONLY result archive remains <package root>/qstab_test_results.zip,
%   containing report.mat and report.txt. Installation/session evidence is
%   report.publicAcceptance inside MAT. An early preflight/runner failure also
%   writes a clearly failed report (zero executed tests is never success).
%   run_qstab_public does not edit startup.m, userpath or saved MATLAB path.
%
%   The pre-existing test_qstab remains available for targeted/external work.
%   This launcher adds no numerical method and changes no solver tolerance.
ip=inputParser;ip.PartialMatching=false;
addParameter(ip,'Examples',true,@localBool);
addParameter(ip,'Figures',false,@localBool);
addParameter(ip,'DependencyScan',true,@localBool);
addParameter(ip,'Verbose',true,@localBool);
parse(ip,varargin{:});o=ip.Results;
root=fileparts(mfilename('fullpath'));
[~,attributes]=fileattrib(root);root=attributes.Name;
writer=@qstab_write_test_report;
oldPath=path;oldFolder=pwd;started=tic;
session=struct('status','NOT_STARTED','pathRestored',false,'folderRestored',false);
try
    [payload,session]=qstab_public_session(root,@()localExecute(root,o));
    report=payload.report;
    report.publicAcceptance=payload.acceptance;
    session.status='RESTORED';
catch ME
    % The session onCleanup restores the caller before this handler runs.
    report=localFailure(root,'PUBLIC_RUN_ERROR',ME);
    report.publicAcceptance=localAcceptance();
    report.publicAcceptance.status='ERROR';
    report.publicAcceptance.message=[ME.identifier ': ' ME.message];
    session.status='EXCEPTION';
end
session.pathRestored=strcmp(path,oldPath);
session.folderRestored=strcmp(pwd,oldFolder);
report.publicAcceptance.session=session;
report.publicAcceptance.pathAndFolderRestored=session.pathRestored&&session.folderRestored;
report.publicAcceptance.totalSeconds=toc(started);
report.publicAcceptance.sourceRunId=report.runId;
if ~report.publicAcceptance.pathAndFolderRestored
    report.publicAcceptance.ok=false;report.publicAcceptance.status='RESTORATION_FAILED';
end
if ~report.publicAcceptance.ok
    report.ok=false;
    if ~ismember(report.status,{'PUBLIC_PREFLIGHT_FAILED','PUBLIC_RUN_ERROR'})
        report.status='PUBLIC_ACCEPTANCE_ISSUES';
    end
end
% test_qstab's first archive is replaced by the same writer/name, now with
% acceptance evidence. No additional archive or third entry is introduced.
report=writer(report,root,fullfile(root,'qstab_test_results'));
fprintf('\nqstab PUBLIC acceptance: %s | regressions %d passed, %d failed, %d skipped\n', ...
    report.publicAcceptance.status,report.passed,report.failed,report.skipped);
fprintf('Caller path and folder restored: %d\n',report.publicAcceptance.pathAndFolderRestored);
if ~isempty(report.publicAcceptance.message),fprintf('%s\n',report.publicAcceptance.message);end
fprintf('SEND THIS ZIP: %s\n',report.archiveFile);
end

function payload=localExecute(root,o)
A=localAcceptance();
A.options=o;
A.before=qstab_check_installation('DependencyScan',o.DependencyScan);
if ~A.before.ok
    report=localFailure(root,'PUBLIC_PREFLIGHT_FAILED',[]);
    report.sourceManifest=A.before.inventory;
    A.status='PREFLIGHT_FAILED';
    A.message='No regression was started. Inspect publicAcceptance.before in report.mat.';
    payload=struct('report',report,'acceptance',A);
    return
end
opts=struct('CheckPrimaryVectors',true,'UseQuaternion',true,'Timing',false);
report=test_qstab('EigenComparison',true,'EigenOptions',opts,'Verbose',o.Verbose, ...
    'DiagnosticExamples',o.Examples,'BallExamples',o.Examples, ...
    'DiscreteExamples',o.Examples,'BoxExamples',o.Examples, ...
    'PaperExamples',o.Examples,'PaperExtensions',o.Examples,'PaperFigures',o.Figures);
A.testSuiteExecuted=true;A.regressionStatus=report.status;
if o.Figures
    try
        A.figureReview=qstab_public_figures(report.paperFigures,root);
    catch ME
        A.figureReview=struct('ok',false,'status','ERROR','message',ME.message,'identifier',ME.identifier);
    end
end
try
    A.after=qstab_check_installation('DependencyScan',false);
    A.sourcesUnchanged=isequaln(A.before.inventory,A.after.inventory);
catch ME
    A.after=struct('ok',false,'status','ERROR','message',ME.message,'identifier',ME.identifier);
end
% A skipped/error optional study need not have legacyRequested. Preserve the
% regression report instead of throwing while inspecting a missing field.
A.noExternalSolverRequested=true;
if isfield(report.eigenComparison,'legacyRequested')
    A.noExternalSolverRequested=~report.eigenComparison.legacyRequested;
end
A.nativeComplete=report.nativeMATLABRun&&report.nativeQuaternionOK&&report.skipped==0;
A.ok=report.ok&&A.before.ok&&A.after.ok&&A.sourcesUnchanged&& ...
    A.noExternalSolverRequested&&A.nativeComplete&&A.figureReview.ok;
if A.ok,A.status='PASSED';
else,A.status='ISSUES';A.message='Inspect regressions and acceptance evidence; nothing was converted to success.';end
payload=struct('report',report,'acceptance',A);
end

function A=localAcceptance()
A=struct('schema','qstab.public-acceptance.1','ok',false,'status','NOT_STARTED', ...
    'options',struct(),'before',struct(),'after',struct(),'testSuiteExecuted',false, ...
    'regressionStatus','','sourcesUnchanged',false,'noExternalSolverRequested',true, ...
    'nativeComplete',false,'pathAndFolderRestored',false,'message','', ...
    'figureReview',struct('status','NOT_REQUESTED','ok',true,'visualInspectionPerformed',false), ...
    'note',['Finite native acceptance in a restricted search path. No claim of ', ...
    'exhaustive dependency detection, validated arithmetic or publication approval.']);
end

function report=localFailure(root,status,ME)
% Match the established report writer's contract even before tests can start.
[~,token]=fileparts(tempname);stamp=datestr(now,30);
env=struct('version',version,'computer',computer,'isOctave',exist('OCTAVE_VERSION','builtin')~=0, ...
    'quaternionPath',which('quaternion'),'qstabPath',fullfile(root,'src','+qstab','isstable.m'));
report=struct('toolbox',struct('name','qstab - Quaternion Stability Toolbox','version','0.18.0'), ...
    'environment',env,'timestamp',stamp,'runId',[stamp '_' token], ...
    'nativeMATLABRun',~env.isOctave,'nativeQuaternionOK',false,'ok',false,'status',status, ...
    'passed',0,'failed',0,'skipped',0,'elapsedSeconds',0,'totalExecutionSeconds',0, ...
    'cases',struct('group',{},'name',{},'status',{},'seconds',{},'message',{},'identifier',{}), ...
    'executionError',struct());
if ~isempty(ME)
    report.executionError=struct('identifier',ME.identifier,'message',ME.message,'stack',ME.stack);
end
end
function tf=localBool(x)
tf=islogical(x)&&isscalar(x);
end
