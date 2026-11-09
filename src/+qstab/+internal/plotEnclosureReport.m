function H = plotEnclosureReport(E,varargin)
%PLOTENCLOSUREREPORT Internal implementation returning the full report.
ip=inputParser;ip.FunctionName='qstab.plotEnclosure';
addParameter(ip,'Parent',[]);addParameter(ip,'Visible','on');
addParameter(ip,'RegionIndices',[]);addParameter(ip,'MaxRegions',128,@localCount);
addParameter(ip,'ShowControls',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Samples',[]);parse(ip,varargin{:});o=ip.Results;
if ~isstruct(E)||~isscalar(E)||~all(isfield(E,{'schema','lambda','regions'}))|| ...
        ~strcmp(E.schema,'qstab.value-enclosure.1')||~iscell(E.regions)||isempty(E.regions)
    error('qstab:plotEnclosure:Input','Supply valueenclosure or certificateEnclosure output.');
end
n=numel(E.regions);idx=o.RegionIndices;if isempty(idx),idx=1:min(n,o.MaxRegions);end
if ~isnumeric(idx)||~isreal(idx)||~isvector(idx)||isempty(idx)||any(~isfinite(idx))|| ...
        any(idx~=fix(idx))||any(idx<1|idx>n)||numel(unique(idx))~=numel(idx)
    error('qstab:plotEnclosure:Indices','RegionIndices must be distinct valid region indices.');
end
for j=idx(:).'
    r=E.regions{j};
    if ~isstruct(r)||~isfield(r,'polygon')||~isnumeric(r.polygon)||isempty(r.polygon)|| ...
            ~isvector(r.polygon)||any(~isfinite(real(r.polygon(:)))|~isfinite(imag(r.polygon(:))))
        error('qstab:plotEnclosure:Polygon','A selected polygon is missing or nonfinite.');
    end
end
samples=[];
if ~isempty(o.Samples)
    S=o.Samples;
    if ~isstruct(S)||~isscalar(S)||~all(isfield(S,{'lambda','points'}))|| ...
            ~isnumeric(S.lambda)||~isnumeric(S.points)||size(S.points,2)~=numel(S.lambda)
        error('qstab:plotEnclosure:Samples','Samples must be valueset output with matching columns.');
    end
    match=find(abs(S.lambda-E.lambda)<=128*eps*max(1,abs(E.lambda)));
    if numel(match)~=1,error('qstab:plotEnclosure:Lambda','Exactly one sample column must match the enclosure argument.');end
    samples=S.points(:,match);
    if any(~isfinite(real(samples))|~isfinite(imag(samples)))
        error('qstab:plotEnclosure:Samples','The selected sample column is nonfinite.');
    end
end
[fig,ax,restore]=qstab.internal.diagnosticAxes(o.Parent,o.Visible,'qstab outer norm-value enclosures'); %#ok<ASGLU>
colors=get(ax,'ColorOrder');outlines=gobjects(numel(idx),1);controlHandles=gobjects(0,1);
for j=1:numel(idx)
    r=E.regions{idx(j)};p=r.polygon(:);color=colors(1+mod(j-1,size(colors,1)),:);
    if numel(p)==1,marker='o';else,marker='none';p=[p;p(1)];end
    outlines(j)=line('Parent',ax,'XData',real(p),'YData',imag(p), ...
        'LineStyle','-','Marker',marker,'LineWidth',1.25,'Color',color, ...
        'DisplayName',sprintf('outer cell %d',idx(j)));
    if o.ShowControls&&isfield(r,'controls')&&~isempty(r.controls)
        c=r.controls(:);controlHandles(end+1,1)=line('Parent',ax,'XData',real(c),'YData',imag(c), ...
            'LineStyle','none','Marker','x','Color',color,'DisplayName','Bernstein controls (not samples)'); %#ok<AGROW>
    end
end
hs=line('Parent',ax,'XData',real(samples),'YData',imag(samples),'LineStyle','none', ...
    'Marker','.','MarkerSize',11,'DisplayName','supplied parameter samples');
hz=line('Parent',ax,'XData',0,'YData',0,'LineStyle','none','Marker','+', ...
    'MarkerSize',12,'LineWidth',1.8,'DisplayName','zero');
handles=outlines(1);labels={'local OUTER sets (union)'};
if ~isempty(controlHandles),handles(end+1)=controlHandles(1);labels{end+1}='controls, not admissible samples';end
if ~isempty(samples),handles(end+1)=hs;labels{end+1}='supplied samples (not a certificate)';end
handles(end+1)=hz;labels{end+1}='zero';legend(ax,handles,labels,'Interpreter','none','Location','best');
xlabel(ax,'Re N_p(lambda)','Interpreter','none');ylabel(ax,'Im N_p(lambda)','Interpreter','none');
grid(ax,'on');axis(ax,'equal');omitted=n-numel(idx);
heading=sprintf('Outer norm-value sets at lambda = %.5g %+.5gi',real(E.lambda),imag(E.lambda));
if isfield(E,'frequencyBounds')&&numel(E.frequencyBounds)==2
    heading=sprintf('Stored whole-frequency outer cell: [%.5g, %.5g]',E.frequencyBounds(1),E.frequencyBounds(2));
end
if isfield(E,'originalVariable')&&strcmp(E.originalVariable,'zinv')
    heading=[heading,' (forward z; input was z^{-1})'];
end
note='Zero in an outer set is inconclusive; this picture is not a stability test.';
if omitted>0,note=sprintf('%s %d regions omitted.',note,omitted);end
title(ax,{heading,note},'Interpreter','none');
H=struct('figure',fig,'axes',ax,'outlines',outlines,'controls',controlHandles,'samples',hs, ...
    'zero',hz,'regionIndices',idx(:).','omittedRegionCount',omitted, ...
    'displaysAllRegions',omitted==0,'samplesMembershipRechecked',false,'isCertificate',false,'data',E);
end
function tf=localCount(x)
tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==fix(x);
end
