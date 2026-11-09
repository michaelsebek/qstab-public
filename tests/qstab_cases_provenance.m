function cases = qstab_cases_provenance(manifestFunction,packageRoot)
%QSTAB_CASES_PROVENANCE Read-only source fingerprints with explicit fallback.
cases={ ...
 'Source fingerprints use SHA-256 when available and never fabricate hashes',@()knownBytes(manifestFunction); ...
 'Source inventory is independent of cwd and excludes runtime results',@()independent(manifestFunction); ...
 'Installed source manifest retains public namespace names and no qeig alias',@()installed(manifestFunction,packageRoot)};
end
function knownBytes(f)
root=tempname;mkdir(root);cleaner=onCleanup(@()rmdir(root,'s')); %#ok<NASGU>
fid=fopen(fullfile(root,'fixture.m'),'wb');assert(fid>=0);fwrite(fid,uint8('abc'),'uint8');fclose(fid);
R=f(root);assert(R.sourceCount==1&&R.sources.bytes==3);
if R.hashingAvailable
 assert(R.complete&&strcmp(R.sources.sha256,'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad'));
else
 assert(~R.complete&&isempty(R.sources.sha256)&&~isempty(R.message));
end
fid=fopen(fullfile(root,'fixture.m'),'rb');raw=fread(fid,Inf,'*uint8');fclose(fid);assert(isequal(raw(:).',uint8('abc')));
end
function independent(f)
root=tempname;mkdir(root);cleaner=onCleanup(@()rmdir(root,'s')); %#ok<NASGU>
mkdir(fullfile(root,'qstab_test_results'));fid=fopen(fullfile(root,'qstab_test_results','ignore.m'),'w');fprintf(fid,'ignored');fclose(fid);
mkdir(fullfile(root,'src'));fid=fopen(fullfile(root,'src','fixture.m'),'w');fprintf(fid,'abc');fclose(fid);
A=f(root);old=pwd;restore=onCleanup(@()cd(old));cd(tempdir); %#ok<NASGU>
B=f(root);assert(isequaln(A,B)&&A.sourceCount==1&&strcmp(A.sources.path,'src/fixture.m'));
end
function installed(f,root)
R=f(root);names={R.sources.path};assert(any(strcmp(names,'src/+qstab/eig.m')));
assert(~any(strcmp(names,'qeig.m'))&&R.sourceCount>=50);
assert(any(strcmp(names,'private/qstab_write_test_report.m')));
end
