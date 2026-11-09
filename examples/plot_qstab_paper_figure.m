function H=plot_qstab_paper_figure(I,varargin)
%PLOT_QSTAB_PAPER_FIGURE Plot one SAVED item; never calls a numerical solver.
% H=plot_qstab_paper_figure(F.items(k),'Visible','on')
% Options: Parent=[] (new figure), Visible='on'/'off', Legend=true.
% Existing axes children and exact NextPlot state are retained. Labels/title
% belong to the selected illustration. No files are created by this function.
qstab_validate_paper_figure(I);
ip=inputParser;ip.FunctionName='plot_qstab_paper_figure';ip.PartialMatching=false;
addParameter(ip,'Parent',[]);addParameter(ip,'Visible','on');
addParameter(ip,'Legend',true,@(x)islogical(x)&&isscalar(x));parse(ip,varargin{:});o=ip.Results;
if isstring(o.Visible)&&isscalar(o.Visible),o.Visible=char(o.Visible);end
if ~ischar(o.Visible)||~ismember(o.Visible,{'on','off'}),error('qstab:paperFigure:Visible','Visible must be on or off.');end
if ~isempty(o.Parent)&&(~isscalar(o.Parent)||~isgraphics(o.Parent,'axes'))
 error('qstab:paperFigure:Parent','Parent must be one axes.');end
created=isempty(o.Parent);
if created
 fig=figure('Visible',o.Visible,'Name',I.id,'NumberTitle','off','Position',[100 100 820 580]);ax=axes('Parent',fig);
else,ax=o.Parent;fig=ancestor(ax,'figure');end
try
 np=get(ax,'NextPlot');cleanup=onCleanup(@()restoreAxes(ax,np)); %#ok<NASGU>
 H=struct('figure',fig,'axes',ax,'lines',{{}},'legend',[],'sourceId',I.id,'isCertificate',false,'data',I);
 palette=get(ax,'ColorOrder');leg={};labels={};
 for j=1:numel(I.curves)
  c=I.curves(j);args={'Parent',ax,'XData',c.x,'YData',c.y, ...
    'LineStyle',c.lineStyle,'Marker',c.marker,'LineWidth',1.35,'MarkerSize',6, ...
    'Color',palette(1+mod(j-1,size(palette,1)),:),'DisplayName',c.label,'Tag',['qstab.paper.',c.role]};
  if I.dimension==3,args=[args,{'ZData',c.z}];end
  h=line(args{:});H.lines{end+1}=h;
  if c.showLegend,leg{end+1}=h;labels{end+1}=c.label;end
 end
 xlabel(ax,I.xlabel,'Interpreter','none');ylabel(ax,I.ylabel,'Interpreter','none');
 if I.dimension==3,zlabel(ax,I.zlabel,'Interpreter','none');view(ax,3);end
 title(ax,[{I.title},wrapCaption(I.subtitle,90)],'Interpreter','none','FontWeight','normal','FontSize',11);
 if I.equalAxes,axis(ax,'equal');end
 set(ax,'Box','on');
 if o.Legend&&~isempty(leg)
  H.legend=legend(ax,[leg{:}],labels,'Interpreter','none','Location','best');
 end
catch ME
 if created&&isgraphics(fig),close(fig);end
 rethrow(ME)
end
end
function restoreAxes(ax,np),if isgraphics(ax,'axes'),set(ax,'NextPlot',np);end,end

function lines=wrapCaption(text,limit)
words=strsplit(text);lines={''};
for j=1:numel(words)
 if isempty(lines{end}),lines{end}=words{j};
 elseif numel(lines{end})+numel(words{j})+1<=limit,lines{end}=[lines{end},' ',words{j}];
 else,lines{end+1}=words{j};end
end
end
