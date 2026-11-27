function report=test_qstab_oneparameter(varargin)
%TEST_QSTAB_ONEPARAMETER Check the new scalar one-parameter functionality.
%   R=test_qstab_oneparameter runs the numerical, article-data and optional
%   exact-rational regression groups. It needs no private toolbox classes.
%   R=test_qstab_oneparameter('ExactCertificates',true) requires usable
%   Symbolic Math Toolbox; false explicitly omits that part; 'auto' (default)
%   records it as pending when Symbolic Math Toolbox is unavailable.
%
%   Options: ExactCertificates='auto', Verbose=true.
%   R.ok means all executed checks passed and explicit requirements were met.
%   R.complete additionally requires all exact checks. Neither field claims
%   native quaternion acceptance; run RUN_QSTAB_PUBLIC for that separate gate.
%   Source examples with native quaternion are run by RUN_QSTAB_ONEPARAMETER.
%   An auxiliary Octave execution is always labelled native-MATLAB pending.
%
%   The current report is saved under <package>/oneparameter_test_results;
%   a new verified <package>/qstab_oneparameter_test_results.zip is made even
%   if a test fails. Its full path is R.archiveFile.
root=fileparts(mfilename('fullpath'));oldPath=path;
restore=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'src'),fullfile(root,'tests'),fullfile(root,'examples'));
ip=inputParser;ip.FunctionName=mfilename;ip.PartialMatching=false;
addParameter(ip,'ExactCertificates','auto',@exactOption);
addParameter(ip,'Verbose',true,@(v)islogical(v)&&isscalar(v));parse(ip,varargin{:});o=ip.Results;
[available,reason]=symbolicAvailable();
requested=~(islogical(o.ExactCertificates)&&~o.ExactCertificates);
required=islogical(o.ExactCertificates)&&o.ExactCertificates;
runExact=requested&&available;
[~,token]=fileparts(tempname);runId=[datestr(now,'yyyymmddTHHMMSS') '_' token];
report=struct('schema','qstab.oneparameter.tests.1','toolbox',qstab.about(), ...
    'runId',runId,'runtime',version,'computer',computer, ...
    'nativeMATLABRun',exist('OCTAVE_VERSION','builtin')==0, ...
    'options',o,'ok',false,'complete',false,'status','RUNNING', ...
    'passed',0,'failed',0,'skipped',0,'elapsedSeconds',0, ...
    'exact',struct('available',available,'executed',runExact,'reason',reason), ...
    'outputDirectory',fullfile(root,'oneparameter_test_results'), ...
    'archiveFile',fullfile(root,'qstab_oneparameter_test_results.zip'), ...
    'cases',struct('group',{},'name',{},'status',{},'seconds',{},'identifier',{},'message',{},'stack',{}));
suites={'numeric',qstab_cases_oneparameter_numeric(); ...
        'examples',qstab_cases_oneparameter_examples(); ...
        'exact',qstab_cases_oneparameter_exact()};
started=tic;
for g=1:size(suites,1)
    cases=suites{g,2};
    for k=1:size(cases,1)
        row=struct('group',suites{g,1},'name',cases{k,1},'status','PASS', ...
            'seconds',0,'identifier','','message','', ...
            'stack',struct('file',{},'name',{},'line',{}));t=tic;
        if strcmp(row.group,'exact')&&~runExact
            row.status='SKIP';row.message=reason;
            if ~requested,row.message='Exact certificates explicitly not requested.';end
            report.skipped=report.skipped+1;
        else
            try,cases{k,2}();report.passed=report.passed+1;
            catch problem
                row.status='FAIL';row.identifier=problem.identifier;row.message=problem.message;
                row.stack=problem.stack;
                report.failed=report.failed+1;
            end
        end
        row.seconds=toc(t);report.cases(end+1)=row;
        if o.Verbose
            fprintf('%s [%s] %s\n',row.status,row.group,row.name);
            if strcmp(row.status,'FAIL')
                fprintf('  %s\n',row.message);
                if ~isempty(row.stack),fprintf('  At %s (line %d).\n',row.stack(1).name,row.stack(1).line);end
            end
        end
    end
end
report.elapsedSeconds=toc(started);
report.ok=report.failed==0&&(~required||available);report.complete=report.ok&&runExact;
if ~report.ok,report.status='REVIEW_REQUIRED';
elseif ~report.complete,report.status='NUMERIC_OK_EXACT_PENDING';
elseif ~report.nativeMATLABRun,report.status='AUXILIARY_OK_NATIVE_PENDING';
else,report.status='ALL_EXECUTED_CHECKS_PASSED';end
writeReport(report);
fprintf('qstab one-parameter tests: %s | %d passed, %d failed, %d skipped\n', ...
    report.status,report.passed,report.failed,report.skipped);
fprintf('SEND THIS ZIP: %s\n',report.archiveFile);
end

function tf=exactOption(v)
tf=(islogical(v)&&isscalar(v))||((ischar(v)||(isstring(v)&&isscalar(v)))&&strcmpi(char(v),'auto'));
end
function [tf,reason]=symbolicAvailable()
tf=false;reason='Symbolic Math Toolbox is unavailable or could not be initialized.';
if exist('sym','file')==2||exist('sym','class')==8
    try,sym(1);tf=true;reason='Symbolic arithmetic is usable.';catch problem,reason=problem.message;end
end
end
function writeReport(report)
out=report.outputDirectory;if ~isfolder(out),mkdir(out);end
save(fullfile(out,'report.mat'),'report','-v7');
fid=fopen(fullfile(out,'report.txt'),'wt');
if fid<0,error('qstab:oneparameter:Report','Cannot write the test report.');end
c=onCleanup(@()fclose(fid));
fprintf(fid,'qstab %s | %s\nRun: %s\nRuntime: %s | %s\n', ...
    report.toolbox.version,report.status,report.runId,report.runtime,report.computer);
fprintf(fid,'Native MATLAB: %d | complete exact scope: %d\n',report.nativeMATLABRun,report.complete);
fprintf(fid,'Passed: %d | failed: %d | skipped: %d | seconds: %.6f\n', ...
    report.passed,report.failed,report.skipped,report.elapsedSeconds);
for k=1:numel(report.cases)
    r=report.cases(k);fprintf(fid,'%s [%s] %s\n',r.status,r.group,r.name);
    if ~isempty(r.message),fprintf(fid,'  %s %s\n',r.identifier,r.message);end
    if strcmp(r.status,'FAIL')
        for j=1:numel(r.stack)
            fprintf(fid,'  %s (line %d): %s\n',r.stack(j).name,r.stack(j).line,r.stack(j).file);
        end
    end
end
fprintf(fid,'Archive: %s\n',report.archiveFile);clear c
staging=tempname;mkdir(staging);cleanup=onCleanup(@()removeStage(staging)); %#ok<NASGU>
z=fullfile(staging,'candidate.zip');entries={'report.mat','report.txt'};
zip(z,entries,out);verify=fullfile(staging,'verify');mkdir(verify);unzip(z,verify);
for k=1:numel(entries)
    if ~isequal(bytes(fullfile(out,entries{k})),bytes(fullfile(verify,entries{k})))
        error('qstab:oneparameter:Archive','Report ZIP read-back did not match the new report.');
    end
end
[ok,msg]=movefile(z,report.archiveFile,'f');
if ~ok,error('qstab:oneparameter:Archive','Cannot replace the result ZIP: %s',msg);end
end
function b=bytes(p)
fid=fopen(p,'rb');if fid<0,error('qstab:oneparameter:Report','Cannot read %s.',p);end
c=onCleanup(@()fclose(fid));b=fread(fid,Inf,'*uint8'); %#ok<NASGU>
end
function removeStage(p)
if isfolder(p),rmdir(p,'s');end
end
