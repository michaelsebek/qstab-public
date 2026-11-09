function H = plotSpectrumReport(lambda, varargin)
%PLOTSPECTRUMREPORT Internal implementation returning the full report.

if nargin < 1
    error('qstab:plotSpectrum:NotEnoughInputs','Supply spectral representatives.');
end
opt = struct('variable','s','dimension',2,'parent',[],'visible','on','tolerance',100*eps);
if mod(numel(varargin),2) ~= 0
    error('qstab:plotSpectrum:NameValue','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = local_text(varargin{k},'OptionName');
    name = lower(name);
    if ~isfield(opt,name)
        error('qstab:plotSpectrum:Option','Unknown option "%s".',name);
    end
    opt.(name) = varargin{k+1};
end
opt.variable = lower(local_text(opt.variable,'Variable'));
if ~any(strcmp(opt.variable,{'s','z','zinv'}))
    error('qstab:plotSpectrum:Variable','Variable must be ''s'', ''z'', or ''zinv''.');
end
if ~isnumeric(opt.dimension) || ~isscalar(opt.dimension) || ...
        ~ismember(opt.dimension,[2 3])
    error('qstab:plotSpectrum:Dimension','Dimension must be 2 or 3.');
end
if ~isnumeric(opt.tolerance) || ~isreal(opt.tolerance) || ...
        ~isscalar(opt.tolerance) || ~isfinite(opt.tolerance) || opt.tolerance < 0
    error('qstab:plotSpectrum:Tolerance','Tolerance must be a finite nonnegative real scalar.');
end
opt.visible = lower(local_text(opt.visible,'Visible'));
if ~any(strcmp(opt.visible,{'on','off'}))
    error('qstab:plotSpectrum:Visible','Visible must be ''on'' or ''off''.');
end
[a,b,c,d] = local_plot_components(lambda);
if ~isequal(size(a),size(b),size(c),size(d)) || ...
        ~isreal(a) || ~isreal(b) || ~isreal(c) || ~isreal(d)
    error('qstab:plotSpectrum:Components','Quaternion components must be equally sized real arrays.');
end
a = double(a(:)); b = double(b(:)); c = double(c(:)); d = double(d(:));
if any(isnan(a)) || any(isnan(b)) || any(isnan(c)) || any(isnan(d))
    error('qstab:plotSpectrum:NaN','NaN representatives cannot be plotted as spectral classes.');
end
rho = hypot(hypot(b,c),d);
finite = isfinite(a) & isfinite(rho);
infinityCount = nnz(~finite);
af = a(finite);
rf = rho(finite);
radius = hypot(af,rf);
if strcmp(opt.variable,'s')
    defect = af;
    tol = opt.tolerance*max(1,radius);
    regionText = 'Stable classes: Re < 0';
elseif strcmp(opt.variable,'z')
    defect = radius-1;
    tol = opt.tolerance*max(1,radius);
    regionText = 'Stable classes: inside the unit boundary';
else
    defect = 1-radius;
    tol = opt.tolerance*max(1,radius);
    regionText = 'Stable classes: outside the unit boundary';
end
status = repmat({'boundary'},numel(af),1);
stable = defect < -tol;
unstable = defect > tol;
status(stable) = {'stable'};
status(unstable) = {'unstable'};
boundary = ~(stable | unstable);

if isempty(opt.parent)
    fig = figure('Color','w','Visible',opt.visible,'Name','qstab spectral classes');
    ax = axes('Parent',fig);
else
    ax = opt.parent;
    if ~isscalar(ax) || ~isgraphics(ax,'axes')
        error('qstab:plotSpectrum:Parent','Parent must be a scalar axes handle.');
    end
    fig = ancestor(ax,'figure');
end
oldNextPlot = get(ax,'NextPlot');
hadContent = ~isempty(get(ax,'Children'));
oldXLim = get(ax,'XLim'); oldYLim = get(ax,'YLim'); oldZLim = get(ax,'ZLim');
restoreHold = onCleanup(@() local_restore_hold(ax,oldNextPlot)); %#ok<NASGU>
set(ax,'NextPlot','add');

[xLimits,yLimits,zLimits] = local_limits(af,rf,opt.variable,opt.dimension);
if hadContent
    xLimits = [min(xLimits(1),oldXLim(1)) max(xLimits(2),oldXLim(2))];
    yLimits = [min(yLimits(1),oldYLim(1)) max(yLimits(2),oldYLim(2))];
    if opt.dimension == 3
        zLimits = [min(zLimits(1),oldZLim(1)) max(zLimits(2),oldZLim(2))];
    end
end
colors = [.12 .53 .28; .80 .13 .16; .85 .54 .08];
markers = {'o','x','d'};
labels = {'stable','unstable','boundary / within tolerance'};
masks = {stable,unstable,boundary};
regionHandle = gobjects(0,1);
classHandles = gobjects(0,1);
legendHandles = gobjects(3,1);

if opt.dimension == 2
    if strcmp(opt.variable,'s')
        if xLimits(1) < 0
            regionHandle = patch(ax,[xLimits(1) 0 0 xLimits(1)], ...
                [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], [.91 .96 .92], ...
                'EdgeColor','none','FaceAlpha',.35,'HandleVisibility','off');
        end
        boundaryHandle = plot(ax,[0 0],yLimits,'--','Color',[.35 .35 .40], ...
            'LineWidth',1.2,'HandleVisibility','off');
    else
        theta = linspace(0,pi,241);
        boundaryHandle = plot(ax,cos(theta),sin(theta),'--','Color',[.35 .35 .40], ...
            'LineWidth',1.2,'HandleVisibility','off');
    end
    for k = 1:3
        % Property-based construction creates exactly one line object even
        % when this stability group has no points (plot may return []).
        legendHandles(k) = line('Parent',ax,'XData',af(masks{k}).', ...
            'YData',rf(masks{k}).','Marker',markers{k}, ...
            'LineStyle','none','Color',colors(k,:),'MarkerSize',7,'LineWidth',1.5, ...
            'DisplayName',labels{k});
    end
    classHandles = legendHandles;
    xlabel(ax,'Re(q)');
    ylabel(ax,'|Im(q)|');
    title(ax,{'Complex representatives of similarity classes',regionText});
    view(ax,2);
else
    if strcmp(opt.variable,'s')
        [yy,zz] = meshgrid(yLimits,zLimits);
        boundaryHandle = surf(ax,zeros(size(yy)),yy,zz, ...
            'FaceColor',[.55 .59 .65],'FaceAlpha',.15, ...
            'EdgeColor',[.50 .53 .58],'LineStyle','--','HandleVisibility','off');
    else
        [xx,yy,zz] = sphere(32);
        boundaryHandle = surf(ax,xx,yy,zz, ...
            'FaceColor',[.55 .65 .80],'FaceAlpha',.10, ...
            'EdgeColor','none','HandleVisibility','off');
        theta = linspace(0,2*pi,241);
        regionHandle = plot3(ax,cos(theta),sin(theta),zeros(size(theta)), ...
            '--','Color',[.35 .40 .50],'HandleVisibility','off');
    end
    theta = linspace(0,2*pi,181);
    for k = 1:3
        which = find(masks{k});
        for j = which(:).'
            if rf(j) == 0
                hj = plot3(ax,af(j),0,0,markers{k},'Color',colors(k,:), ...
                    'MarkerSize',7,'LineWidth',1.5,'HandleVisibility','off');
            else
                hj = plot3(ax,af(j)*ones(size(theta)),rf(j)*cos(theta),rf(j)*sin(theta), ...
                    'Color',colors(k,:),'LineWidth',1.5,'HandleVisibility','off');
            end
            classHandles(end+1,1) = hj; %#ok<AGROW>
        end
        % Keep a scalar legend handle also for an empty stability group.
        legendHandles(k) = line('Parent',ax,'XData',af(masks{k}).', ...
            'YData',rf(masks{k}).','ZData',zeros(1,nnz(masks{k})), ...
            'Marker',markers{k},'LineStyle','none','Color',colors(k,:),'MarkerSize',6, ...
            'LineWidth',1.5,'DisplayName',labels{k});
    end
    xlabel(ax,'Re(q)');
    ylabel(ax,'Im_i(q)');
    zlabel(ax,'Im_j(q)');
    title(ax,{'Exact class slice Im_k(q) = 0',regionText});
    view(ax,38,24);
    zlim(ax,zLimits);
end
xlim(ax,xLimits);
ylim(ax,yLimits);
daspect(ax,[1 1 1]);
grid(ax,'on');
box(ax,'on');
legend(ax,legendHandles,labels,'Location','best');
annotationHandle = gobjects(0,1);
if infinityCount > 0
    annotationHandle = text(ax,'Units','normalized','Position',[.02 .98 0], ...
        'String',sprintf('Infinite representatives: %d (shown separately)',infinityCount), ...
        'VerticalAlignment','top','Interpreter','none');
end

H = struct();
H.figure = fig;
H.axes = ax;
H.classes = classHandles;
H.markers = legendHandles;
H.boundary = boundaryHandle;
H.region = regionHandle;
H.annotation = annotationHandle;
H.variable = opt.variable;
H.dimension = opt.dimension;
H.canonicalLambda = complex(af,rf);
H.originalIndices = find(finite);
H.classification = status;
H.stableMask = stable;
H.unstableMask = unstable;
H.boundaryMask = boundary;
H.infinityCount = infinityCount;
H.infinityIndices = find(~finite);
H.tolerance = opt.tolerance;
H.isCertificate = false;
H.geometry = 'exact class quotient';
if opt.dimension == 3, H.geometry = 'exact slice Im_k(q)=0'; end
end

function s = local_text(s,kind)
if isstring(s) && isscalar(s), s = char(s); end
if ~ischar(s) || size(s,1) ~= 1
    error(['qstab:plotSpectrum:' kind],'Text options must be character rows or scalar strings.');
end
end

function [xl,yl,zl] = local_limits(a,r,variable,dimension)
if strcmp(variable,'s')
    landmarks = [0; a];
else
    landmarks = [-1; 1; a];
end
lo = min(landmarks); hi = max(landmarks);
span = max(1,hi-lo);
xl = [lo-.12*span hi+.12*span];
rmax = 1.12*max([1; r]);
if dimension == 2
    yl = [0 rmax];
else
    yl = [-rmax rmax];
end
zl = [-rmax rmax];
end

function local_restore_hold(ax,value)
if isgraphics(ax,'axes'), set(ax,'NextPlot',value); end
end

function [a,b,c,d] = local_plot_components(Q)
% Core qstab.components intentionally rejects Inf. Plotting needs to retain
% infinite roots for a separate annotation, so only this reader allows Inf.
if isa(Q,'quaternion')
    [a,b,c,d] = parts(Q);
elseif isnumeric(Q) || islogical(Q)
    a = real(Q); b = imag(Q); c = zeros(size(Q)); d = zeros(size(Q));
elseif isstruct(Q) && isscalar(Q) && all(isfield(Q,{'a','b','c','d'}))
    a = Q.a; b = Q.b; c = Q.c; d = Q.d;
    if ~isnumeric(a) || ~isnumeric(b) || ~isnumeric(c) || ~isnumeric(d)
        error('qstab:plotSpectrum:Components','Component fields must be numeric.');
    end
else
    error('qstab:plotSpectrum:Input', ...
        'Use numeric, native quaternion, or component-struct representatives.');
end
end
