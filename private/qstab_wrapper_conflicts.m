function conflicts=qstab_wrapper_conflicts(root)
%QSTAB_WRAPPER_CONFLICTS Read-only check before enabling unqualified aliases.
%   Only names actually distributed in this package's qwrappers are checked.
%   A same-named foreign function anywhere on the current MATLAB path is a
%   conflict, even when the qstab alias currently precedes it. No function is
%   executed and no file/path is changed. Missing/off-path qfiles cannot be
%   audited by this check; a source inventory is a separate operation.
%
%   qeig is deliberately NOT distributed. qstab.eig remains namespaced.

    wrapperDir=fullfile(root,'qwrappers');
    entries=dir(fullfile(wrapperDir,'*.m'));
    conflicts=struct('name',{},'locations',{});
    for k=1:numel(entries)
        [~,name]=fileparts(entries(k).name);
        own=fullfile(wrapperDir,entries(k).name);
        hits=which(name,'-all');
        if isempty(hits), continue; end
        if ischar(hits), hits=cellstr(hits); end
        if isstring(hits), hits=cellstr(hits); end
        foreign={};
        for j=1:numel(hits)
            hit=char(hits{j});
            if ~isempty(hit) && ~localSameFile(hit,own)
                foreign{end+1,1}=hit; %#ok<AGROW>
            end
        end
        if ~isempty(foreign)
            conflicts(end+1)=struct('name',name,'locations',{foreign}); %#ok<AGROW>
        end
    end
end

function same=localSameFile(a,b)
    a=strrep(a,'\','/'); b=strrep(b,'\','/');
    if ispc, same=strcmpi(a,b); else, same=strcmp(a,b); end
end
