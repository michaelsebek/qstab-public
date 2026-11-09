function [D,v,g,h,unbounded] = familyPlotModel(F,o)
%FAMILYPLOTMODEL Validate public inputs before sampling complete members.
kind=o.kind; v=o.variable; g=o.grade; h=[]; unbounded=false;
if any(strcmp(kind,{'one-parameter','oneparam'})),kind='oneparameter';end
tagged=isstruct(F)&&isscalar(F)&&isfield(F,'schema');
if tagged && isequal(F.schema,'qstab.family-definition.1')
    [raw,k,fv,fg]=qstab.internal.familyDefinition(F);
    kind=agree(kind,k,'Uncertainty'); v=agree(v,fv,'Variable'); g=agree(g,fg,'Grade');
    h=F.samplingPeriod;
elseif tagged && isequal(F.schema,'qstab.scalar-family.1')
    % Reconstruct a public definition, then validate that dependent fields
    % describe that same model. Never silently accept stale edited caches.
    [raw,k]=restoreRaw(F); kind=agree(kind,k,'Uncertainty');
    if isfield(F,'variable') && ~isempty(F.variable),v=agree(v,F.variable,'Variable');end
    if isfield(F,'samplingPeriod'),h=F.samplingPeriod;end
    if isfield(F,'grade')&&~isempty(F.grade),g=agree(g,F.grade,'Grade');end
elseif tagged
    error('qstab:sample:Family','Unknown family schema; see help qstab.familycomponents.');
else
    raw=F;
end
if isstruct(raw)&&isscalar(raw)&&all(isfield(raw,{'p0','p1'})) && ...
        (~isfield(raw,'bounds') || isequal(raw.bounds,[]))
    if isempty(o.box)
        error('qstab:sample:SampleBounds','An unrestricted quaternion has no finite uniform sampling domain. Supply ParameterBox as four finite [lower upper] rows.');
    end
    unbounded=true; raw.bounds=o.box;
end
D=qstab.familycomponents(raw,'Uncertainty',kind);
if tagged && isequal(F.schema,'qstab.scalar-family.1')
    keys={'kind','degree','base','generators','bounds','vertices','side','intervalLower','intervalUpper','parameterIndex'};
    for n=1:numel(keys)
        k=keys{n};
        if ~isfield(F,k)||~isequal(F.(k),D.(k))
            error('qstab:sample:Family','Inconsistent derived field %s. Rebuild the family from its definition.',k);
        end
    end
end
if isempty(v),v='s';end
if isempty(g),g=D.degree;end
if ~(ischar(v)&&any(strcmp(v,{'s','z','zinv'})))
    error('qstab:sample:Variable','Variable must be s, z or zinv.');
end
if ~(isnumeric(g)&&isscalar(g)&&isreal(g)&&isfinite(g)&&g==fix(g)&&g>=D.degree)
    error('qstab:sample:Grade','Grade cannot be smaller than the supplied coefficient layout.');
end
end

function value=agree(request,stored,key)
if ~isempty(request)&&~isequal(request,stored)
    error('qstab:sample:Family','Explicit %s conflicts with stored family metadata.',key);
end
value=stored;
end

function [raw,kind]=restoreRaw(F)
keys={'degree','kind','base','generators','vertices','bounds','intervalLower','intervalUpper','side'};
if ~all(isfield(F,keys)) || ~isnumeric(F.degree)||~isscalar(F.degree)||~isreal(F.degree)|| ...
    ~isfinite(F.degree)||F.degree<0||F.degree~=floor(F.degree)
    error('qstab:sample:Family','Use a complete qstab.familycomponents result or pass the original family.');
end
kind=F.kind;
switch kind
    case 'polytope'
        V=F.vertices;
        if ~isnumeric(V)||~isreal(V)||ndims(V)>3||isempty(V)||size(V,1)~=4||size(V,2)~=F.degree+1
            error('qstab:sample:Family','Invalid canonical polytope vertices.');
        end
        raw=cell(1,size(V,3));for j=1:numel(raw),raw{j}=pack(V(:,:,j),F.degree);end
    case 'interval'
        raw={pack(F.intervalLower,F.degree),pack(F.intervalUpper,F.degree)};
    case {'affine','oneparameter'}
        G=F.generators;
        if ~isnumeric(G)||~isreal(G)||ndims(G)>3||size(G,1)~=4||size(G,2)~=F.degree+1
            error('qstab:sample:Family','Invalid canonical affine generators.');
        end
        if strcmp(kind,'oneparameter')
            if size(G,3)~=4,error('qstab:sample:Family','Four generators are required for a quaternion parameter.');end
            raw=struct('p0',pack(F.base,F.degree),'p1',pack(G(:,:,1),F.degree),'side',F.side,'bounds',F.bounds);
        else
            generators=cell(1,size(G,3));
            for j=1:numel(generators),generators{j}=pack(G(:,:,j),F.degree);end
            raw=struct('center',pack(F.base,F.degree),'generators',{generators},'bounds',F.bounds);
        end
    otherwise
        error('qstab:sample:Family','Unknown canonical family kind.');
end
end

function p=pack(T,g)
if ~isnumeric(T)||~isreal(T)||~isequal(size(T),[4,g+1])||any(~isfinite(T(:)))
    error('qstab:sample:Family','Each coefficient page must be a finite real 4-by-(grade+1) array.');
end
p=struct('a',T(1,:),'b',T(2,:),'c',T(3,:),'d',T(4,:));
end
