function H = plotValueSetReport(S, varargin)
%PLOTVALUESETREPORT Internal implementation returning the full report.

if ~isstruct(S) || ~isscalar(S) || ~isfield(S,'points') || ...
        ~isnumeric(S.points) || ndims(S.points) > 2 || isempty(S.points)
    error('qstab:plotValueSet:Input','S must contain a nonempty numeric matrix S.points.');
end
opt = struct('parent',[],'visible','on','lambdaindex',1:size(S.points,2),'connect',false);
if mod(numel(varargin),2) ~= 0
    error('qstab:plotValueSet:NameValue','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = local_option_name(varargin{k});
    if ~isfield(opt,name)
        error('qstab:plotValueSet:Option','Unknown option "%s".',name);
    end
    opt.(name) = varargin{k+1};
end
opt.visible = local_visible(opt.visible);
idx = opt.lambdaindex;
if ~isnumeric(idx) || ~isreal(idx) || ~isvector(idx) || isempty(idx) || ...
        any(~isfinite(idx(:))) || any(idx(:) ~= fix(idx(:))) || ...
        any(idx(:) < 1) || any(idx(:) > size(S.points,2))
    error('qstab:plotValueSet:LambdaIndex','LambdaIndex must contain valid positive column indices.');
end
if ~isscalar(opt.connect) || ~(islogical(opt.connect) || isnumeric(opt.connect)) || ...
        ~isreal(opt.connect) || ~isfinite(opt.connect) || ~ismember(opt.connect,[0 1])
    error('qstab:plotValueSet:Connect','Connect must be true or false.');
end
if isempty(opt.parent)
    fig = figure('Color','w','Visible',opt.visible,'Name','qstab norm-polynomial values');
    ax = axes('Parent',fig);
else
    ax = opt.parent;
    if ~isscalar(ax) || ~isgraphics(ax,'axes')
        error('qstab:plotValueSet:Parent','Parent must be a scalar axes handle.');
    end
    fig = ancestor(ax,'figure');
end
oldNextPlot = get(ax,'NextPlot');
restoreHold = onCleanup(@() local_restore_hold(ax,oldNextPlot)); %#ok<NASGU>
set(ax,'NextPlot','add');

points = double(S.points(:,idx));
finite = isfinite(real(points)) & isfinite(imag(points));
points(~finite) = complex(NaN,NaN);
if logical(opt.connect)
    style = '.-';
else
    style = '.';
end
if size(points,2) > 8 && ~logical(opt.connect)
    hp = plot(ax,real(points(:)),imag(points(:)),style,'MarkerSize',13, ...
        'LineWidth',1,'DisplayName','sampled values');
else
    hp = gobjects(size(points,2),1);
    for k = 1:size(points,2)
        hp(k) = plot(ax,real(points(:,k)),imag(points(:,k)),style, ...
            'MarkerSize',13,'LineWidth',1, ...
            'DisplayName',sprintf('samples, column %d',idx(k)));
    end
end
hz = plot(ax,0,0,'+','Color',[0.80 0.12 0.14],'MarkerSize',12,'LineWidth',2, ...
    'DisplayName','zero');
xlabel(ax,'Re N_p(\lambda)');
ylabel(ax,'Im N_p(\lambda)');
grid(ax,'on');
axis(ax,'equal');
if isfield(S,'kind') && strcmp(S.kind,'singleton')
    title(ax,{'Norm-polynomial values of a fixed polynomial','Values at the supplied complex arguments'});
else
    title(ax,{'Sampled norm-polynomial values','Samples do not certify zero exclusion'});
end
if numel(hp) <= 8
    legend(ax,[hp(:);hz],'Location','best');
else
    legend(ax,[hp(1);hz],{'sampled values','zero'},'Location','best');
end
hnote = gobjects(0,1);
omitted = nnz(~finite);
if omitted > 0
    hnote = text(ax,'Units','normalized','Position',[.02 .98 0], ...
        'String',sprintf('%d nonfinite values omitted',omitted), ...
        'VerticalAlignment','top','Interpreter','none');
end
H = struct('figure',fig,'axes',ax,'points',hp,'zero',hz,'annotation',hnote, ...
    'lambdaIndex',idx(:).','omittedNonfiniteCount',omitted,'isCertificate',false);
end

function name = local_option_name(x)
if isstring(x) && isscalar(x), x = char(x); end
if ~ischar(x) || size(x,1) ~= 1
    error('qstab:plotValueSet:OptionName','Option names must be character rows or scalar strings.');
end
name = lower(x);
end

function v = local_visible(v)
if isstring(v) && isscalar(v), v = char(v); end
if ~ischar(v) || ~any(strcmpi(v,{'on','off'}))
    error('qstab:plotValueSet:Visible','Visible must be ''on'' or ''off''.');
end
v = lower(v);
end

function local_restore_hold(ax,value)
if isgraphics(ax,'axes'), set(ax,'NextPlot',value); end
end
