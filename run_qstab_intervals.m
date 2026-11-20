function report=run_qstab_intervals(varargin)
%RUN_QSTAB_INTERVALS Reproduce the interval article and its supplement.
%   R = run_qstab_intervals runs all 24 numbered numerical example groups.
%   R = run_qstab_intervals('I19') runs only the five-gain regulator.
%   R = run_qstab_intervals({'I09','I10'}) runs both quartic contacts.
%   C = run_qstab_intervals('list') returns the catalogue and prints its keys.
%   R = run_qstab_intervals('all','Plots',true) also saves eight figure panels.
%
%   Name/value options:
%     UseQuaternion true   Use native MATLAB quaternion inputs. False uses
%                          ordinary component structs, with the same core.
%     Plots         false  Export computed results to PNG/PDF/FIG and CSV.
%     Visible       'off'  'on' displays new figure windows as well.
%     Exact         false  Also recompute exact rational certificates using
%                          Symbolic Math Toolbox. Optional; may take longer.
%     Save          true   Save MAT, text, CSV and a newly created report ZIP.
%     Verbose       true   Print one summary per example and result locations.
%     OutputDir     []     Raw result parent; default <package>/interval_results.
%   Every saved run uses its own subdirectory. The archive is always
%   <package>/qstab_interval_results.zip and contains only this run.
%   Its absolute name is R.archiveFile. A calculation error inside an example
%   is recorded and archived. Invalid options or a missing requested native
%   quaternion class stop at preflight, before a result directory is created.
%
%   Sources are IJRCN_manuscript_v7 and IJRCN_Zenodo-supplement_v5; their
%   exact file hashes and coverage map are in docs/interval_article.md.
%   Calculated results and paper references are separate. R.ok means that
%   the selected numerical comparisons passed; it does not silently turn
%   floating-point diagnostics or plotted samples into exact certificates.
%   R.exact records whether an exact proof was freshly executed in this run.
%
%   Example for the normal PUBLIC environment:
%     qstab_setup('Examples',true);
%     R = run_qstab_intervals('all','Plots',true);
%     assert(R.ok,'Inspect %s',R.archiveFile);
%     R.cases(19).output.radius
%
%   The caller's MATLAB path, current directory, RNG, and existing figures
%   are preserved. No private quaternion toolbox is needed or initialized.
%   See also test_qstab, run_qstab_public, qstab_interval_cases.
root=fileparts(mfilename('fullpath'));oldpath=path;restore=onCleanup(@()path(oldpath)); %#ok<NASGU>
addpath(fullfile(root,'src'),fullfile(root,'examples'));
C=qstab_interval_cases();selection={C.id};
if ~isempty(varargin)
 first=varargin{1};
 if (ischar(first)&&any(strcmpi(first,{'list','all'})))||iscell(first)||(ischar(first)&&numel(first)==3&&first(1)=='I')
  varargin(1)=[];
  if ischar(first)&&strcmpi(first,'list')
   if ~isempty(varargin),error('qstab:intervalStudy:Options','The list command takes no options.');end
   report=C;for k=1:numel(C),fprintf('%s  %s\n',C(k).id,C(k).title);end;return
  elseif ischar(first)&&strcmpi(first,'all')
  elseif ischar(first),selection={first};else,selection=first;end
 end
end
if ~iscell(selection)||isempty(selection)||~all(cellfun(@ischar,selection))||any(~ismember(selection,{C.id}))||numel(unique(selection))~=numel(selection)
 error('qstab:intervalStudy:Selection','Use unique catalogue keys, for example I03 or {''I09'',''I10''}.');
end
C=C(ismember({C.id},selection));
ip=inputParser;ip.FunctionName='run_qstab_intervals';ip.PartialMatching=false;
addParameter(ip,'UseQuaternion',true,@flag);addParameter(ip,'Plots',false,@flag);
addParameter(ip,'Exact',false,@flag);addParameter(ip,'Save',true,@flag);addParameter(ip,'Verbose',true,@flag);
addParameter(ip,'Visible','off',@(x)ischar(x)&&any(strcmp(x,{'on','off'})));
addParameter(ip,'OutputDir','',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
parse(ip,varargin{:});o=ip.Results;
if o.Plots&&~o.Save,error('qstab:intervalStudy:Options','Plots requires Save=true so the figures have a documented destination.');end
if o.UseQuaternion
 try,q=quaternion(0,0,0,0);[a,b,c,d]=parts(q); %#ok<ASGLU>
 catch err,error('qstab:intervalStudy:Quaternion','Native quaternion unavailable. Use UseQuaternion=false for a component-only run. %s',err.message);end
end
[~,token]=fileparts(tempname);runId=[datestr(now,30) '_' token];
out=char(o.OutputDir);if isempty(out),out=fullfile(root,'interval_results');
elseif ~absolute(out),out=fullfile(root,out);end
out=fullfile(out,runId);
report=struct('schema','qstab.interval-study.1','toolbox',qstab.about(), ...
 'runId',runId,'runtime',version,'computer',computer,'isOctave',exist('OCTAVE_VERSION','builtin')~=0, ...
 'nativeQuaternionInputs',o.UseQuaternion,'options',o,'sourceCatalogue',C, ...
 'ok',true,'status','RUNNING','passed',0,'failed',0,'checksPassed',0,'checksFailed',0, ...
 'elapsedSeconds',0,'outputDirectory',out,'archiveFile',fullfile(root,'qstab_interval_results.zip'), ...
 'cases',struct('id',{},'title',{},'sourceAnchor',{},'status',{},'output',{},'checks',{},'seconds',{},'error',{}), ...
 'exact',struct('executed',false,'status','NOT_REQUESTED','ok',false), ...
 'graphics',struct('executed',false,'files',{{}},'error',''), ...
 'assurance','Numerical comparisons are distinct from exact continuum proof replay.');
started=tic;
for k=1:numel(C)
 t=tic;row=struct('id',C(k).id,'title',C(k).title,'sourceAnchor',C(k).sourceAnchor, ...
  'status','ERROR','output',struct(),'checks',struct([]),'seconds',0,'error',struct());
 try
  [row.output,row.checks]=qstab_interval_evaluate(row.id,o);
  if ~isempty(row.checks)&&all([row.checks.pass]),row.status='MATCHED';report.passed=report.passed+1;
  else,row.status='MISMATCH';report.failed=report.failed+1;end
  report.checksPassed=report.checksPassed+sum([row.checks.pass]);report.checksFailed=report.checksFailed+sum(~[row.checks.pass]);
 catch err
  row.error=struct('identifier',err.identifier,'message',err.message,'stack',err.stack);report.failed=report.failed+1;
 end
 row.seconds=toc(t);report.cases(end+1)=row;
 if o.Verbose
  fprintf('%s %-8s %s (%.3f s)\n',row.id,row.status,row.title,row.seconds);
  if strcmp(row.status,'ERROR'),fprintf('  %s\n',row.error.message);
  elseif strcmp(row.status,'MISMATCH'),for j=find(~[row.checks.pass]),fprintf('  MISMATCH: %s\n',row.checks(j).name);end;end
 end
end
if o.Exact
 try,report.exact=verify_qstab_interval_certificates('Verbose',o.Verbose);report.exact.executed=true;
 catch err,report.exact=struct('executed',true,'ok',false,'status','ERROR','identifier',err.identifier,'message',err.message);end
end
if o.Save
 if ~isfolder(out),mkdir(out);end
 if o.Plots
  try,report.graphics.files=plot_qstab_intervals(report,'Directory',out,'Visible',o.Visible);report.graphics.executed=true;
  catch err,report.graphics.executed=true;report.graphics.error=err.message;end
 end
end
report.elapsedSeconds=toc(started);report.ok=report.failed==0&&isempty(report.graphics.error)&&(~o.Exact||report.exact.ok);
if report.ok,report.status='NUMERICAL_EXAMPLES_MATCHED';else,report.status='REVIEW_REQUIRED';end
if o.Exact&&report.ok,report.status='NUMERICAL_AND_EXACT_CHECKS_PASSED';end
if o.Save
 save(fullfile(out,'report.mat'),'report','-v7');writeText(report,fullfile(out,'report.txt'));writeCSV(report,fullfile(out,'comparison.csv'));
 entries={'report.mat','report.txt','comparison.csv'};
 for k=1:numel(report.graphics.files),[~,name,ext]=fileparts(report.graphics.files{k});entries{end+1}=[name ext];end %#ok<AGROW>
 if isfile(report.archiveFile),delete(report.archiveFile);end
 zip(report.archiveFile,entries,out);
 if ~isfile(report.archiveFile),error('qstab:intervalStudy:Archive','The result archive was not created.');end
end
if o.Verbose
 fprintf('\nqstab interval article: %s | %d groups matched, %d failed | %d checks passed\n',report.status,report.passed,report.failed,report.checksPassed);
 fprintf('Exact rational replay: %s\n',report.exact.status);
 if o.Save,fprintf('Results: %s\nSEND THIS ZIP: %s\n',out,report.archiveFile);end
end
end
function y=flag(x),y=islogical(x)&&isscalar(x);end
function y=absolute(x),y=~isempty(x)&&(x(1)=='/'||x(1)=='\'||(numel(x)>1&&x(2)==':'));end
function writeText(R,file)
f=fopen(file,'w');if f<0,error('qstab:intervalStudy:Write','Cannot write %s',file);end;c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'qstab interval article and supplement\nVersion: %s\nRun: %s\nRuntime: %s | %s\nStatus: %s\nNative quaternion inputs: %d\n',R.toolbox.version,R.runId,R.runtime,R.computer,R.status,R.nativeQuaternionInputs);
fprintf(f,'Exact rational replay: %s\n%s\n\n',R.exact.status,R.assurance);
for k=1:numel(R.cases)
 A=R.cases(k);fprintf(f,'%s %s | %s\nSource: %s\n',A.id,A.status,A.title,A.sourceAnchor);
 for j=1:numel(A.checks),B=A.checks(j);fprintf(f,'  pass=%d %s\n    computed: %s\n    reference: %s\n',B.pass,B.name,mat2str(B.actual,16),mat2str(B.reference,16));end
 if strcmp(A.status,'ERROR'),fprintf(f,'  ERROR %s: %s\n',A.error.identifier,A.error.message);end
 fprintf(f,'\n');
end
fprintf(f,'Elapsed seconds: %.6f\nArchive: %s\n',R.elapsedSeconds,R.archiveFile);
end
function writeCSV(R,file)
f=fopen(file,'w');if f<0,error('qstab:intervalStudy:Write','Cannot write %s',file);end;c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'case,check,passed,computed,reference,tolerance\n');
for k=1:numel(R.cases),A=R.cases(k);for j=1:numel(A.checks),B=A.checks(j);fprintf(f,'%s,"%s",%d,"%s","%s",%.17g\n',A.id,strrep(B.name,'"','""'),B.pass,mat2str(B.actual,16),mat2str(B.reference,16),B.tolerance);end;end
end
