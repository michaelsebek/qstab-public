function [n,o]=familyPlotOptions(args,valuesMode)
% Parse the two plotting interfaces. N is a model count, never per dimension.
n=100;
if ~isempty(args) && (isnumeric(args{1}) || islogical(args{1})), n=args{1}; args(1)=[]; end
if ~(isnumeric(n)&&isscalar(n)&&isreal(n)&&isfinite(n)&&n>=1&&n==fix(n)&&n<=1e6)
    error('qstab:sample:Count','N must be an integer in [1,1e6]: the number of sampled polynomials.');
end
n=double(n);
o=struct('kind','','variable','','grade',[],'sampling','grid','seed',[],'box',[],'plot',true,'parent',[], ...
    'visible','on','axis',[1 0 0],'sphere',false,'directions',32,'projection','class');
seen={}; hasK=false;
if mod(numel(args),2)~=0,error('qstab:sample:Options','Use complete name/value pairs after the model count.');end
for h=1:2:numel(args)
    t=textOption(args{h}); v=args{h+1};
    if any(strcmp(t,seen)), error('qstab:sample:Options','Option "%s" is duplicated.',t); end
    seen{end+1}=t; %#ok<AGROW>
    switch t
        case 'uncertainty', o.kind=textOption(v);
        case 'variable'
            o.variable=textOption(v);
            if any(strcmp(o.variable,{'z^-1','delay'})), o.variable='zinv'; end
            if ~any(strcmp(o.variable,{'s','z','zinv'})), error('qstab:sample:Variable','Variable must be s, z or zinv.'); end
        case 'grade'
            if ~(isnumeric(v)&&isreal(v)&&isscalar(v)&&isfinite(v)&&v>=0&&v==fix(v)),error('qstab:sample:Grade','Grade must be a nonnegative integer.');end
            o.grade=double(v);
        case 'sampling', o.sampling=textOption(v);
        case 'seed', o.seed=v;
        case 'parameterbox', o.box=v;
        case 'plot'
            if ~(islogical(v)&&isscalar(v)), error('qstab:sample:Options','Plot must be true or false.');end
            o.plot=v;
        case 'visible'
            o.visible=textOption(v);
            if ~any(strcmp(o.visible,{'on','off'})),error('qstab:sample:Options','Visible must be on or off.');end
        case 'parent'
            if ~isempty(v) && ~(isscalar(v)&&isgraphics(v,'axes')),error('qstab:sample:Options','Parent must be a scalar axes handle.');end
            o.parent=v;
        case 'axisunit'
            if ~valuesMode,error('qstab:sample:Options','AxisUnit applies to plotFamilyValues only.');end
            if (ischar(v) || isstring(v)) && strcmp(textOption(v),'sphere')
                o.sphere=true;
            else
                o.axis=unit(v);
            end
        case 'directioncount'
            if ~valuesMode,error('qstab:sample:Options','DirectionCount applies to plotFamilyValues only.');end
            if ~(isnumeric(v)&&isscalar(v)&&isreal(v)&&isfinite(v)&&v==fix(v)&&v>=1&&v<=1e6)
                error('qstab:sample:Directions','DirectionCount must be an integer from 1 to 1e6.');
            end
            hasK=true;o.directions=double(v);
        case 'projection'
            if ~valuesMode,error('qstab:sample:Options','Projection applies to plotFamilyValues only.');end
            o.projection=textOption(v);
            if ~any(strcmp(o.projection,{'class','slice','normpoly'})),error('qstab:sample:Projection','Use class, slice or normpoly.');end
        otherwise,error('qstab:sample:Options','Unknown option "%s".',t);
    end
end
if o.sphere
    if ~strcmp(o.projection,'class')
        error('qstab:sample:Projection','Use sphere with class. For slice specify one signed axis; for normpoly omit sphere (no orientation search is needed).');
    end
    if n>floor(1e6/o.directions)
        error('qstab:sample:Directions','Reduce N or DirectionCount: at most 1e6 member/direction evaluations per sphere plot.');
    end
    o.axis=sphericalGrid(o.directions);
elseif hasK
    error('qstab:sample:Directions','DirectionCount applies only when AxisUnit is sphere.');
end
end
function t=textOption(v)
if isstring(v)&&isscalar(v),v=char(v);end
if ~(ischar(v)&&isrow(v)&&~isempty(v)),error('qstab:sample:Options','Expected scalar text for an option.');end
t=lower(strtrim(v));
end
function u=unit(v)
if ischar(v)||isstring(v)
    t=textOption(v); j=find(strcmp(t,{'i','j','k'}),1);
    if isempty(j),error('qstab:sample:Axis','AxisUnit must be i, j, k, a pure unit quaternion, a real unit 3-vector or sphere.');end
    u=zeros(1,3);u(j)=1;return
elseif isnumeric(v)&&isreal(v)&&isvector(v)&&numel(v)==3
    u=double(v(:).');
else
    if ~(isa(v,'quaternion')||(isnumeric(v)&&isscalar(v))||(isstruct(v)&&isscalar(v)&&all(isfield(v,{'a','b','c','d'}))))
        error('qstab:sample:Axis','AxisUnit must be one pure unit quaternion or a real unit 3-vector.');
    end
    [a,b,c,d]=qstab.components(v);
    if ~isscalar(a)||a~=0,error('qstab:sample:Axis','AxisUnit requires a single quaternion with zero real part.');end
    u=[b c d];
end
l=hypot(hypot(u(1),u(2)),u(3));
if ~isfinite(l)||abs(l-1)>1e-7,error('qstab:sample:Axis','AxisUnit must have unit length.');end
u=u/l;
end
function U=sphericalGrid(k)
anchors=[eye(3);-eye(3)];
if k<=6,U=anchors(1:k,:);return;end
m=k-6; j=(0:m-1).'; h=1-2*(j+.5)/m;
phi=(pi*(3-sqrt(5)))*(j+.5);
r=sqrt(max(0,1-h.^2)); U=[anchors; r.*cos(phi), r.*sin(phi), h];
end
