function R = qstab_source_manifest(root)
%QSTAB_SOURCE_MANIFEST Read-only fingerprints of installed MATLAB sources.
% Diagnostic metadata only: an unavailable SHA-256 backend does not prevent
% test execution or create a fabricated digest. No third report/ZIP entry is
% introduced. Runtime results folders and version-control metadata are ignored.
R=struct('algorithm','SHA-256','hashingAvailable',false,'complete',false, ...
    'sourceCount',0,'sources',struct('path',{},'sha256',{},'bytes',{},'error',{}), ...
    'message','');
paths=localWalk(root,'');paths=sort(paths);R.sourceCount=numel(paths);
try
    md=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
    R.hashingAvailable=true;
catch ME
    md=[];R.message=['SHA-256 backend unavailable: ',ME.message];
end
for k=1:numel(paths)
    rel=paths{k};entry=struct('path',strrep(rel,'\','/'),'sha256','', ...
        'bytes',NaN,'error','');
    fid=fopen(fullfile(root,rel),'rb');
    if fid<0
        entry.error='Source file could not be read.';
    else
        closer=onCleanup(@()fclose(fid));
        raw=fread(fid,Inf,'*uint8');entry.bytes=numel(raw);
        clear closer
        if R.hashingAvailable
            try
                md.reset();md.update(typecast(raw(:),'int8'));
                digest=typecast(md.digest(),'uint8');
                entry.sha256=lower(reshape(dec2hex(digest,2).',1,[]));
            catch ME
                entry.error=ME.message;
            end
        else
            entry.error=R.message;
        end
    end
    R.sources(end+1)=entry; %#ok<AGROW>
end
R.complete=R.hashingAvailable&&all(cellfun(@isempty,{R.sources.error}))&& ...
    all(cellfun(@(x)numel(x)==64,{R.sources.sha256}));
if R.complete,R.message='All listed MATLAB source files were fingerprinted without modification.';
elseif isempty(R.message),R.message='Some source files could not be fingerprinted.';end
end
function files=localWalk(root,relative)
files={};entries=dir(fullfile(root,relative));
for j=1:numel(entries)
    name=entries(j).name;
    if ismember(name,{'.','..','.git','.svn','qstab_test_results'}),continue;end
    rel=fullfile(relative,name);
    if entries(j).isdir
        files=[files,localWalk(root,rel)]; %#ok<AGROW>
    else
        [~,~,ext]=fileparts(name);
        if strcmpi(ext,'.m'),files{end+1}=rel;end %#ok<AGROW>
    end
end
end
