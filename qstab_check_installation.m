function A=qstab_check_installation(varargin)
%QSTAB_CHECK_INSTALLATION Read-only integrity and MATLAB resolution checks.
%   A=qstab_check_installation verifies all installed .m files against the
%   distributed MANIFEST.sha256 and checks public/example function paths.
%   Missing, modified, unexpected and unhashed code are reported separately.
%   It never changes MATLAB path, repairs files, initializes a toolbox or
%   calls a solver. In particular it never calls an external qeig.
%
%   Options (exact names):
%     PackageRoot      this file's directory by default
%     CheckResolution  true; require selected functions in this installation
%                      and the native quaternion constructor under matlabroot
%     DependencyScan   false; optionally call MATLAB's static dependency tool
%                      for src/+qstab public entries (not the external study)
%
%   A.ok is an installation check, not a numerical or publication certificate.
%   A.integrity lists source problems; A.resolution lists actual paths.
%   Static dependencies cannot prove completeness of dynamic dispatch and
%   list only detected installed products. A failed/unavailable requested scan
%   is not silently converted to success. SHA-256 is not a signed manifest.
%   For temporary isolation and execution use RUN_QSTAB_PUBLIC.
%
%   CheckResolution=false is a file-only inspection of another package root;
%   no source file in that other tree is executed. DependencyScan requires
%   CheckResolution=true. All outputs are ordinary, serializable data.
ip=inputParser;ip.PartialMatching=false;
addParameter(ip,'PackageRoot',fileparts(mfilename('fullpath')),@localText);
addParameter(ip,'CheckResolution',true,@localBool);
addParameter(ip,'DependencyScan',false,@localBool);
parse(ip,varargin{:});o=ip.Results;
if o.DependencyScan&&~o.CheckResolution
    error('qstab:PublicAuditOptions','DependencyScan requires CheckResolution=true.');
end
[found,rootInfo]=fileattrib(char(o.PackageRoot));
if ~found||~isscalar(rootInfo)||~rootInfo.directory
    error('qstab:PublicAuditRoot','PackageRoot must be an existing directory.');
end
root=rootInfo.Name;M=qstab_source_manifest(root);
manifest=fullfile(root,'MANIFEST.sha256');text='';readError='';
try,text=fileread(manifest);catch ME,readError=ME.message;end
I=qstab_public_manifest(text,M);
if ~isempty(readError),I.errors{end+1}=readError;I.ok=false;end
A=struct('schema','qstab.installation-check.1','packageRoot',root,'ok',false, ...
    'integrity',I,'inventory',M,'resolution',localEmptyResolution(), ...
    'dependencyScan',localEmptyDependencies(),'matlabVersion',version, ...
    'matlabRoot',matlabroot,'computer',computer,'note', ...
    'Finite installation/dependency checks, not proof of every possible execution branch.');
if o.CheckResolution
    names={'qstab_setup','test_qstab','run_qstab_public','qstab_check_installation'};
    expected=cellfun(@(n)fullfile(root,[n '.m']),names,'UniformOutput',false);
    for k=1:numel(I.expected)
        rel=I.expected(k).path;t=regexp(rel,'^src/\+qstab/([^/]+)\.m$','tokens','once');
        if ~isempty(t)&&~strcmpi(t{1},'Contents')
            names{end+1}=['qstab.',t{1}];expected{end+1}=fullfile(root,rel); %#ok<AGROW>
        elseif strncmp(rel,'examples/',9)&&numel(strfind(rel,'/'))==1
            [~,base]=fileparts(rel);
            names{end+1}=base;expected{end+1}=fullfile(root,rel); %#ok<AGROW>
        end
    end
    rows=struct('name',{},'expected',{},'selected',{},'allLocations',{},'ok',{});
    for k=1:numel(names)
        selected=which(names{k});locations=which(names{k},'-all');
        locations=localLocations(locations);
        same=localSame(selected,expected{k});
        uniqueOwn=all(cellfun(@(p)localSame(p,expected{k}),locations));
        rows(end+1)=struct('name',names{k},'expected',expected{k}, ...
            'selected',selected,'allLocations',{locations},'ok',same&&uniqueOwn); %#ok<AGROW>
    end
    qpath=which('quaternion');qOK=localWithin(qpath,matlabroot);
    A.resolution=struct('status','CHECKED','ok',all([rows.ok])&&qOK, ...
        'entries',rows,'quaternionPath',qpath,'nativeQuaternionLocation',qOK, ...
        'legacyVisible',which('qeig'),'legacyRequested',false);
end
if o.DependencyScan
    if I.ok&&A.resolution.ok
        entryFiles={};
        for k=1:numel(I.expected)
            rel=I.expected(k).path;
            if ~isempty(regexp(rel,'^src/\+qstab/[^/]+\.m$','once'))&& ...
                    isempty(regexp(rel,'/Contents\.m$','once'))
                entryFiles{end+1}=fullfile(root,rel); %#ok<AGROW>
            end
        end
        try
            [fList,pList]=matlab.codetools.requiredFilesAndProducts(entryFiles);
            fList=localLocations(fList);external={};
            for k=1:numel(fList)
                if ~(localWithin(fList{k},root)||localWithin(fList{k},matlabroot))
                    external{end+1}=fList{k}; %#ok<AGROW>
                end
            end
            A.dependencyScan=struct('status','COMPLETED','ok',isempty(external), ...
                'entryFiles',{entryFiles},'files',{fList},'products',pList, ...
                'externalFiles',{external},'message', ...
                'Static detected dependencies only; dynamic calls and uninstalled products can be missed.');
        catch ME
            A.dependencyScan=localEmptyDependencies();
            A.dependencyScan.status='ERROR';A.dependencyScan.ok=false;
            A.dependencyScan.message=[ME.identifier ': ' ME.message];
        end
    else
        A.dependencyScan.status='SKIPPED_INSTALLATION_ISSUES';A.dependencyScan.ok=false;
    end
end
A.ok=I.ok&&(~o.CheckResolution||A.resolution.ok)&& ...
    (~o.DependencyScan||A.dependencyScan.ok);
if nargout==0
    fprintf('qstab installation: ok=%d | expected %d, installed %d MATLAB sources\n', ...
        A.ok,I.expectedCount,I.observedCount);
    fprintf('Missing %d | modified %d | unexpected %d | unverified %d\n', ...
        numel(I.missing),numel(I.modified),numel(I.unexpected),numel(I.unverified));
    fprintf('Resolution: %s | dependencies: %s\n',A.resolution.status,A.dependencyScan.status);
end
end
function R=localEmptyResolution()
R=struct('status','NOT_REQUESTED','ok',false,'entries', ...
    struct('name',{},'expected',{},'selected',{},'allLocations',{},'ok',{}), ...
    'quaternionPath','','nativeQuaternionLocation',false,'legacyVisible','', ...
    'legacyRequested',false);
end
function D=localEmptyDependencies()
D=struct('status','NOT_REQUESTED','ok',true,'entryFiles',{{}},'files',{{}}, ...
    'products',struct([]),'externalFiles',{{}},'message','');
end
function c=localLocations(x)
if isempty(x),c={};elseif iscell(x),c=x(:).';elseif isstring(x),c=cellstr(x(:)).';
else,c=cellstr(x);c=c(:).';end
end
function tf=localSame(a,b)
a=strrep(a,'\','/');b=strrep(b,'\','/');
if ispc,tf=strcmpi(a,b);else,tf=strcmp(a,b);end
end
function tf=localWithin(p,root)
p=strrep(p,'\','/');root=strrep(root,'\','/');
while numel(root)>1&&root(end)=='/',root(end)=[];end
if ispc,p=lower(p);root=lower(root);end
tf=~isempty(p)&&(strcmp(p,root)||strncmp(p,[root '/'],numel(root)+1));
end
function tf=localText(x)
tf=(ischar(x)&&isrow(x)&&~isempty(strtrim(x)))|| ...
    (isstring(x)&&isscalar(x)&&~ismissing(x)&&strlength(x)>0);
end
function tf=localBool(x)
tf=islogical(x)&&isscalar(x);
end
