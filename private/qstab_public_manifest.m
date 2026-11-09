function R=qstab_public_manifest(text,inventory,ignoreCase)
%QSTAB_PUBLIC_MANIFEST Compare installed MATLAB sources with a release manifest.
% Pure data check: never opens paths named by the manifest, executes code,
% repairs files, or changes path. SHA-256 consistency is not authentication.
if nargin<3,ignoreCase=ispc;end
if ~(ischar(text)&&isstruct(inventory)&&isscalar(inventory)&& ...
        isfield(inventory,'sources')&&isfield(inventory,'complete'))
    error('qstab:PublicManifestInput','Text and a source inventory are required.');
end
R=struct('schema','qstab.public-manifest.1','ok',false,'manifestValid',false, ...
    'inventoryComplete',logical(inventory.complete),'expectedCount',0, ...
    'observedCount',numel(inventory.sources),'missing',{{}},'unexpected',{{}}, ...
    'modified',{{}},'unverified',{{}},'errors',{{}},'expected', ...
    struct('path',{},'sha256',{}),'note', ...
    'Unsigned SHA-256 manifest checks file consistency, not publisher identity.');
lines=regexp(text,'\r\n|\n|\r','split');
seen={};expectedKeys={};
for k=1:numel(lines)
    line=strtrim(lines{k});
    if isempty(line)||line(1)=='#',continue;end
    t=regexp(line,'^([0-9A-Fa-f]{64})[ \t]+\*?(.+)$','tokens','once');
    if isempty(t)
        R.errors{end+1}=sprintf('Malformed manifest line %d.',k); %#ok<AGROW>
        continue
    end
    rel=strrep(strtrim(t{2}),'\','/');
    % Preserve empty components: src//file.m must not become a valid path.
    pieces=strsplit(rel,'/','CollapseDelimiters',false);
    if isempty(rel)||rel(1)=='/'||any(rel==':')|| ...
            any(cellfun(@isempty,pieces))||any(ismember(pieces,{'.','..'}))
        R.errors{end+1}=sprintf('Unsafe manifest path on line %d.',k); %#ok<AGROW>
        continue
    end
    key=localKey(rel,ignoreCase);
    if any(strcmp(seen,key))
        R.errors{end+1}=['Duplicate manifest path: ',rel]; %#ok<AGROW>
        continue
    end
    seen{end+1}=key; %#ok<AGROW>
    [~,~,ext]=fileparts(rel);
    if strcmpi(ext,'.m')
        R.expected(end+1)=struct('path',rel,'sha256',lower(t{1})); %#ok<AGROW>
        expectedKeys{end+1}=key; %#ok<AGROW>
    end
end
R.expectedCount=numel(R.expected);
if R.expectedCount==0,R.errors{end+1}='Manifest contains no MATLAB sources.';end
R.manifestValid=isempty(R.errors);
actual=inventory.sources;actualKeys=cell(1,numel(actual));
for k=1:numel(actual)
    actualKeys{k}=localKey(actual(k).path,ignoreCase);
end
if numel(unique(actualKeys))~=numel(actualKeys)
    R.errors{end+1}='Installed inventory contains duplicate source paths.';
end
for k=1:numel(R.expected)
    ix=find(strcmp(actualKeys,expectedKeys{k}));
    if isempty(ix)
        R.missing{end+1}=R.expected(k).path; %#ok<AGROW>
    elseif numel(ix)==1
        value=actual(ix).sha256;
        if ~isempty(actual(ix).error)||isempty(regexp(value,'^[0-9a-fA-F]{64}$','once'))
            R.unverified{end+1}=R.expected(k).path; %#ok<AGROW>
        elseif ~strcmpi(value,R.expected(k).sha256)
            R.modified{end+1}=R.expected(k).path; %#ok<AGROW>
        end
    end
end
for k=1:numel(actual)
    if ~any(strcmp(expectedKeys,actualKeys{k}))
        R.unexpected{end+1}=actual(k).path; %#ok<AGROW>
    end
end
R.ok=R.manifestValid&&R.inventoryComplete&&isempty(R.errors)&& ...
    isempty(R.missing)&&isempty(R.unexpected)&&isempty(R.modified)&&isempty(R.unverified);
end
function key=localKey(value,ignoreCase)
key=strrep(value,'\','/');
if ignoreCase,key=lower(key);end
end
