function cases=qstab_cases_public_artifacts(capture,writer)
%QSTAB_CASES_PUBLIC_ARTIFACTS Byte transport, NOT image-decoder validation.
cases={ ...
 'PUBLIC figure review copies original PNG/PDF bytes without rendering',@()copyBytes(capture); ...
 'PUBLIC figure review refuses missing, changed or wrong-signature exports',@()badFiles(capture); ...
 'PUBLIC figure review cannot read arbitrary files through metadata paths',@()badOrigin(capture); ...
 'PNG/PDF byte payloads round-trip through the same two-entry result archive',@()roundtrip(capture,writer)};
end
function [root,F,guard]=fixture()
root=tempname;mkdir(root);mkdir(fullfile(root,'qstab_paper_figures'));
guard=onCleanup(@()rmdir(root,'s'));
% Signature-bearing opaque byte fixtures; no claim of valid image content.
bytes={uint8([137 80 78 71 13 10 26 10 1 2 3]),uint8('%PDF-opaque fixture')};formats={'png','pdf'};
files=struct('id',{},'format',{},'file',{},'bytes',{});
for k=1:2
 p=fullfile(root,'qstab_paper_figures',['figure_a.' formats{k}]);
 fid=fopen(p,'wb');assert(fid>=0);fwrite(fid,bytes{k},'uint8');fclose(fid);
 files(end+1)=struct('id','figure_a','format',formats{k},'file',p,'bytes',numel(bytes{k})); %#ok<AGROW>
end
F=struct('ok',true,'items',struct('id','figure_a'),'export',struct('ok',true,'files',files));
end
function copyBytes(f)
[r,F,g]=fixture(); %#ok<ASGLU,NASGU>
p=path;d=pwd;C=f(F,r);assert(C.ok&&~C.visualInspectionPerformed&&numel(C.files)==2);
for k=1:2
 fid=fopen(F.export.files(k).file,'rb');b=fread(fid,Inf,'*uint8');fclose(fid);
 assert(isequal(C.files(k).data,b)&&isa(C.files(k).data,'uint8'));
end
assert(C.totalBytes==sum([F.export.files.bytes])&&strcmp(path,p)&&strcmp(pwd,d));
end
function badFiles(f)
[r,F,g]=fixture(); %#ok<ASGLU,NASGU>
B=F;B.export.files(1).bytes=B.export.files(1).bytes+1;mustReject(@()f(B,r));
fid=fopen(F.export.files(1).file,'r+b');assert(fid>=0);fwrite(fid,0,'uint8');fclose(fid);
mustReject(@()f(F,r));delete(F.export.files(1).file);mustReject(@()f(F,r));
B=F;B.ok=false;mustReject(@()f(B,r));
end
function badOrigin(f)
[r,F,g]=fixture(); %#ok<ASGLU,NASGU>
B=F;B.export.files(1).file=fullfile(r,'other.png');mustReject(@()f(B,r));
B=F;B.items.id='../other';mustReject(@()f(B,r));
B=F;B.export.files(3)=B.export.files(1);mustReject(@()f(B,r));
end
function roundtrip(f,writer)
[r,F,g]=fixture(); %#ok<ASGLU,NASGU>
C=f(F,r);
R=struct('toolbox',struct('version','fixture'),'environment', ...
 struct('version','fixture','computer','fixture','quaternionPath','','qstabPath',''), ...
 'timestamp','fixture','runId','figure-fixture','nativeMATLABRun',false, ...
 'nativeQuaternionOK',false,'ok',true,'status','FIXTURE','passed',1,'failed',0,'skipped',0,'elapsedSeconds',0, ...
 'cases',struct('group',{},'name',{},'status',{},'seconds',{},'message',{},'identifier',{}), ...
 'publicAcceptance',struct('figureReview',C));
R=writer(R,r,'results');out=fullfile(r,'unpacked');mkdir(out);unzip(R.archiveFile,out);
B=load(fullfile(out,'report.mat'),'report');assert(isequaln(B.report.publicAcceptance.figureReview,C));
list=dir(out);list=list(~ismember({list.name},{'.','..'}));assert(numel(list)==2);
assert(exist(fullfile(out,'report.mat'),'file')==2&&exist(fullfile(out,'report.txt'),'file')==2);
end
function mustReject(f)
failed=false;try,f();catch ME,failed=strcmp(ME.identifier,'qstab:PublicFigures');end
assert(failed,'Unsafe or inconsistent export payload must be rejected.');
end
