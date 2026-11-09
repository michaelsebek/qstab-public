function S=qstab_public_path_plan(currentPath,root,mathWorksRoot,ignoreCase)
%QSTAB_PUBLIC_PATH_PLAN Plan temporary MATLAB-root-only path isolation.
% Pure string policy. No restoredefaultpath/userpath/savepath or filesystem
% mutation. Only the package root/src/examples and existing matlabroot paths
% survive. Product/add-on paths outside matlabroot are deliberately excluded.
if nargin<4,ignoreCase=ispc;end
if ischar(currentPath),dirs=strsplit(currentPath,pathsep);
elseif iscell(currentPath)&&all(cellfun(@ischar,currentPath)),dirs=currentPath;
else,error('qstab:PublicPathInput','Expected a path string or cell array.');end
if ~localAbsolute(root)||~localAbsolute(mathWorksRoot)
    error('qstab:PublicPathInput','Package and MathWorks roots must be absolute.');
end
pkg={root,fullfile(root,'src'),fullfile(root,'examples')};
keep={};removed={};keys={};
for k=1:numel(dirs)
    item=dirs{k};if isempty(item),continue;end
    if localWithin(item,mathWorksRoot,ignoreCase)
        key=localKey(item,ignoreCase);
        if ~any(strcmp(keys,key)),keep{end+1}=item;keys{end+1}=key;end %#ok<AGROW>
    elseif ~any(cellfun(@(p)strcmp(localKey(item,ignoreCase),localKey(p,ignoreCase)),pkg))
        removed{end+1}=item; %#ok<AGROW>
    end
end
S=struct('schema','qstab.public-path-plan.1','packageRoot',root, ...
    'mathWorksRoot',mathWorksRoot,'packagePaths',{pkg},'mathWorksPaths',{keep}, ...
    'removedPaths',{removed},'isolatedPaths',{[pkg,keep]}, ...
    'note','Temporary path restriction, not a process or security sandbox.');
end
function tf=localAbsolute(p)
if ~(ischar(p)&&isrow(p)&&~isempty(p)),tf=false;return;end
q=strrep(p,'\','/');
tf=q(1)=='/'||~isempty(regexp(q,'^[A-Za-z]:/','once'));
end
function tf=localWithin(p,root,ignoreCase)
if ~localAbsolute(p),tf=false;return;end
parts=strsplit(strrep(p,'\','/'),'/');
if any(ismember(parts,{'.','..'})),tf=false;return;end
p=localKey(p,ignoreCase);root=localKey(root,ignoreCase);
tf=strcmp(p,root)||strncmp(p,[root,'/'],numel(root)+1);
end
function k=localKey(p,ignoreCase)
k=strrep(p,'\','/');
while numel(k)>1&&k(end)=='/',k(end)=[];end
if ignoreCase,k=lower(k);end
end
