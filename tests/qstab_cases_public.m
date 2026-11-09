function cases=qstab_cases_public(checkManifest,pathPlan,session,writer,root)
%QSTAB_CASES_PUBLIC Installation consistency and reversible PUBLIC execution.
% No solver is replaced; fixtures are ordinary data or temporary directories.
cases={ ...
 'A complete source inventory matches its unsigned SHA-256 manifest',@()complete(checkManifest); ...
 'Missing example code is reported before test execution',@()missing(checkManifest); ...
 'Modified MATLAB code is distinguished from a missing file',@()modified(checkManifest); ...
 'A stale docs/examples source is an unexpected duplicate',@()unexpected(checkManifest); ...
 'Unavailable or incomplete hashes cannot attest source integrity',@()unhashed(checkManifest); ...
 'Malformed manifest digests are rejected without reading their paths',@()badDigest(checkManifest); ...
 'An absent or source-free manifest is not a successful empty audit',@()emptyManifest(checkManifest); ...
 'Parent and current-directory traversal are invalid manifest paths',@()traversal(checkManifest); ...
 'Absolute drive and UNC manifest paths are rejected',@()absoluteManifest(checkManifest); ...
 'Repeated manifest paths are rejected rather than silently overwritten',@()duplicate(checkManifest); ...
 'Duplicate installed-source records cannot hide an integrity mismatch',@()duplicateInventory(checkManifest); ...
 'Manifest slash and case rules are explicit for Windows and POSIX',@()caseRules(checkManifest); ...
 'PUBLIC path keeps only this root/src/examples and existing MATLAB folders',@()pathSelection(pathPlan); ...
 'A same-prefix sibling directory is not inside matlabroot',@()prefixSafety(pathPlan); ...
 'Development, wrapper and duplicate package directories are not initial PUBLIC paths',@()packagePaths(pathPlan); ...
 'Windows path comparison handles drive case and trailing separators',@()windowsPaths(pathPlan); ...
 'Relative and dot-traversal paths are not accepted as MATLAB folders',@()relativePaths(pathPlan); ...
 'MATLAB path order is preserved and duplicate MATLAB entries are collapsed',@()pathOrder(pathPlan); ...
 'A successful isolated session restores the exact caller path and directory',@()sessionSuccess(session); ...
 'An exception in an isolated session restores the exact caller environment',@()sessionFailure(session); ...
 'Path isolation does not run a startup script in the inspected directory',@()noStartup(session); ...
 'Foreign current-folder functions are invisible during isolated execution',@()foreignCwd(session); ...
 'PUBLIC launcher rejects external-solver overrides before changing paths',@invalidOptions; ...
 'Installed read-only audit survives the unchanged two-entry report ZIP',@()installedRoundtrip(writer,root)};
end

function [text,I]=fixture()
h=repmat('a',1,64);g=repmat('b',1,64);
text=sprintf('%s  src/+qstab/eig.m\n%s  examples/demo.m\n',h,g);
s=struct('path','src/+qstab/eig.m','sha256',h,'bytes',3,'error','');
s(2)=struct('path','examples/demo.m','sha256',g,'bytes',4,'error','');
I=struct('complete',true,'sources',s);
end
function complete(f)
[t,I]=fixture();R=f(t,I,false);assert(R.ok&&R.expectedCount==2&&R.observedCount==2);
assert(isempty(R.errors)&&isempty(R.modified)&&isempty(R.missing));
end
function missing(f)
[t,I]=fixture();I.sources(2)=[];R=f(t,I,false);
assert(~R.ok&&isequal(R.missing,{'examples/demo.m'})&&isempty(R.modified));
end
function modified(f)
[t,I]=fixture();I.sources(2).sha256=repmat('c',1,64);R=f(t,I,false);
assert(~R.ok&&isequal(R.modified,{'examples/demo.m'})&&isempty(R.missing));
end
function unexpected(f)
[t,I]=fixture();I.sources(3)=I.sources(2);I.sources(3).path='docs/examples/demo.m';
R=f(t,I,false);assert(~R.ok&&isequal(R.unexpected,{'docs/examples/demo.m'}));
end
function unhashed(f)
[t,I]=fixture();I.sources(1).sha256='';I.sources(1).error='SHA backend unavailable';I.complete=false;
R=f(t,I,false);assert(~R.ok&&~R.inventoryComplete&&numel(R.unverified)==1&&isempty(R.modified));
[t,I]=fixture();I.complete=false;R=f(t,I,false);assert(~R.ok);
end
function badDigest(f)
[t,I]=fixture();t(1)='g';R=f(t,I,false);assert(~R.ok&&~R.manifestValid&&~isempty(R.errors));
end
function emptyManifest(f)
[~,I]=fixture();R=f('',I,false);assert(~R.ok&&~R.manifestValid);
R=f([repmat('a',1,64) '  README.md'],I,false);assert(~R.ok&&R.expectedCount==0);
end
function traversal(f)
[t,I]=fixture();bad={'../qeig.m','src/../qeig.m','./qeig.m','src//evil.m'};
for k=1:numel(bad)
 R=f([t sprintf('%s  %s\n',repmat('a',1,64),bad{k})],I,false);assert(~R.ok&&~R.manifestValid);
end
end
function absoluteManifest(f)
[t,I]=fixture();bad={'/tmp/qeig.m','C:/tmp/qeig.m','\\server\share\qeig.m','C:qeig.m'};
for k=1:numel(bad)
 R=f([t sprintf('%s  %s\n',repmat('a',1,64),bad{k})],I,false);assert(~R.ok&&~R.manifestValid);
end
end
function duplicate(f)
[t,I]=fixture();R=f([t t],I,false);assert(~R.ok&&~R.manifestValid);
end
function duplicateInventory(f)
[t,I]=fixture();I.sources(3)=I.sources(1);R=f(t,I,false);assert(~R.ok&&~isempty(R.errors));
end
function caseRules(f)
[t,I]=fixture();I.sources(1).path='SRC\+QSTAB\EIG.M';
R=f(t,I,true);assert(R.ok);
R=f(t,I,false);assert(~R.ok&&numel(R.missing)==1&&numel(R.unexpected)==1);
end
function pathSelection(f)
S=f({'/mw/toolbox/base','/third/qfiles','/pkg','/pkg/src','/pkg/examples'},'/pkg','/mw',false);
assert(isequal(S.isolatedPaths,{'/pkg',fullfile('/pkg','src'),fullfile('/pkg','examples'),'/mw/toolbox/base'}));
assert(isequal(S.removedPaths,{'/third/qfiles'}));
end
function prefixSafety(f)
S=f({'/mw_elsewhere/toolbox','/mw2','/mw/toolbox'},'/pkg','/mw',false);
assert(numel(S.mathWorksPaths)==1&&numel(S.removedPaths)==2);
end
function packagePaths(f)
S=f({'/pkg/qwrappers','/pkg/tests','/pkg/dev','/pkg/docs/examples','/other/src'},'/pkg','/mw',false);
assert(numel(S.isolatedPaths)==3&&numel(S.removedPaths)==5);
end
function windowsPaths(f)
S=f({'c:\MATLAB\toolbox\','C:\MATLAB_OLD\toolbox','c:\QSTAB\'},'C:\qstab','C:\MATLAB',true);
assert(numel(S.mathWorksPaths)==1&&numel(S.removedPaths)==1);
end
function relativePaths(f)
S=f({'./mw/toolbox','/mw/../private','/mw/./toolbox','relative','/mw/good'},'/pkg','/mw',false);
assert(isequal(S.mathWorksPaths,{'/mw/good'})&&numel(S.removedPaths)==4);
end
function pathOrder(f)
S=f({'/mw/b','/mw/a','/mw/b','/mw/a/'},'/pkg','/mw',false);
assert(isequal(S.mathWorksPaths,{'/mw/b','/mw/a'}));
end
function [r,g]=temporaryRoot()
r=tempname;mkdir(r);mkdir(fullfile(r,'src'));mkdir(fullfile(r,'examples'));
[~,a]=fileattrib(r);r=a.Name;
g=onCleanup(@()rmdir(r,'s'));
end
function sessionSuccess(f)
[r,g]=temporaryRoot(); %#ok<ASGLU,NASGU>
p=path;d=pwd;[v,S]=f(r,@localObserve);
assert(S.pathRestored&&S.folderRestored&&strcmp(path,p)&&strcmp(pwd,d));
assert(strcmp(v.folder,r)&&any(strcmp(v.paths,r)));
assert(~any(strcmp(v.paths,fullfile(r,'tests'))));
end
function v=localObserve()
v=struct('folder',pwd,'paths',{strsplit(path,pathsep)},'qeig',which('qeig'));
end
function sessionFailure(f)
[r,g]=temporaryRoot(); %#ok<ASGLU,NASGU>
p=path;d=pwd;thrown=false;
try,f(r,@localThrow);catch ME,thrown=strcmp(ME.identifier,'qstab:public:Fixture');end
assert(thrown&&strcmp(path,p)&&strcmp(pwd,d));
end
function v=localThrow() %#ok<STOUT>
error('qstab:public:Fixture','Intentional isolated callback exception.');
end
function noStartup(f)
[r,g]=temporaryRoot(); %#ok<ASGLU,NASGU>
fid=fopen(fullfile(r,'startup.m'),'wt');assert(fid>=0);
fprintf(fid,'error(''qstab:public:StartupRan'',''Must not initialize startup.'');\n');fclose(fid);
[v,S]=f(r,@localObserve);assert(S.pathRestored&&strcmp(v.folder,r));
end
function foreignCwd(f)
[r,g]=temporaryRoot(); %#ok<ASGLU,NASGU>
foreign=tempname;mkdir(foreign);old=pwd;
restore=onCleanup(@()localRestoreFolder(old,foreign)); %#ok<NASGU>
fid=fopen(fullfile(foreign,'qeig.m'),'wt');assert(fid>=0);
fprintf(fid,'function z=qeig(varargin)\nerror(''qstab:public:ForeignCalled'',''Not permitted.'');\nend\n');fclose(fid);
cd(foreign);foreignFolder=pwd;[v,S]=f(r,@localObserve);
assert(isempty(v.qeig)&&S.folderRestored&&strcmp(pwd,foreignFolder));
end
function localRestoreFolder(old,folder)
cd(old);rmdir(folder,'s');
end
function invalidOptions()
p=path;d=pwd;bad={{'LegacyFunction',@sin},{'EigenOptions',struct()},{'Examples',1},{'DependencyScan',[]}};
for k=1:numel(bad)
 rejected=false;try,run_qstab_public(bad{k}{:});catch,rejected=true;end
 assert(rejected&&strcmp(path,p)&&strcmp(pwd,d));
end
end
function installedRoundtrip(writer,root)
p=path;d=pwd;A=qstab_check_installation('PackageRoot',root,'CheckResolution',false);
assert(strcmp(path,p)&&strcmp(pwd,d));
assert(A.integrity.expectedCount==A.inventory.sourceCount);
if A.inventory.hashingAvailable,assert(A.ok);else,assert(~A.ok);end
[r,g]=temporaryRoot(); %#ok<ASGLU,NASGU>
R=struct('toolbox',struct('version','fixture'),'environment', ...
 struct('version','fixture','computer','fixture','quaternionPath','','qstabPath',''), ...
 'timestamp','fixture','runId','public-fixture','nativeMATLABRun',false, ...
 'nativeQuaternionOK',false,'ok',false,'status','PUBLIC_PREFLIGHT_FAILED', ...
 'passed',0,'failed',0,'skipped',0,'elapsedSeconds',0, ...
 'cases',struct('group',{},'name',{},'status',{},'seconds',{},'message',{},'identifier',{}), ...
 'publicAcceptance',struct('ok',false,'before',A));
R=writer(R,r,'results');tmp=fullfile(r,'unpacked');mkdir(tmp);unzip(R.archiveFile,tmp);
B=load(fullfile(tmp,'report.mat'),'report');assert(isequaln(B.report.publicAcceptance,R.publicAcceptance));
contents=dir(tmp);contents=contents(~ismember({contents.name},{'.','..'}));
assert(isequal(sort({contents.name}),{'report.mat','report.txt'}));
end
