function H = plotCertificateReport(R,varargin)
%PLOTCERTIFICATEREPORT Internal implementation returning the full report.
ip=inputParser;ip.FunctionName='qstab.plotCertificate';
addParameter(ip,'View','auto',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'Parent',[]);addParameter(ip,'Visible','on');
addParameter(ip,'RowIndices',[]);addParameter(ip,'MaxRows',200,@localCount);
parse(ip,varargin{:});o=ip.Results;view=lower(char(o.View));
if ~ismember(view,{'auto','coverage','separation'}),error('qstab:plotCertificate:View','View must be auto, coverage, or separation.');end
if isstruct(R)&&isscalar(R)&&isfield(R,'schema')&&strcmp(R.schema,'qstab.certificate-data.1'),D=R;
else,D=qstab.internal.certificateDataReport(R);end
n=numel(D.rows);idx=o.RowIndices;if isempty(idx),idx=1:min(n,o.MaxRows);end
if ~isnumeric(idx)||~isreal(idx)||(~isempty(idx)&&~isvector(idx))||any(~isfinite(idx))|| ...
        any(idx~=fix(idx))||any(idx<1|idx>n)||numel(unique(idx))~=numel(idx)
    error('qstab:plotCertificate:Indices','RowIndices must select distinct stored rows.');
end
if strcmp(view,'auto')
    if n>0&&any(arrayfun(@(r)~isempty(r.frequencyBounds),D.rows)),view='coverage';else,view='separation';end
end
units={};
if strcmp(view,'separation')
    for j=idx(:).'
        if isfinite(D.rows(j).metric),units{end+1}=D.rows(j).metricUnits;end %#ok<AGROW>
    end
    units=unique(units);
    if numel(units)>1,error('qstab:plotCertificate:Units','Metrics with different units must not share one axis.');end
end
[fig,ax,restore]=qstab.internal.diagnosticAxes(o.Parent,o.Visible,'qstab retained certificate diagnostics'); %#ok<ASGLU>
accepted=gobjects(0,1);unresolved=gobjects(0,1);noCoordinate=0;
for k=1:numel(idx)
    j=idx(k);row=D.rows(j);
    if strcmp(view,'coverage')
        if isempty(row.frequencyBounds),noCoordinate=noCoordinate+1;continue;end
        x=row.frequencyBounds;y=[j j];
    else
        if ~isfinite(row.metric),noCoordinate=noCoordinate+1;continue;end
        x=j;y=row.metric;
    end
    if strcmp(row.state,'accepted'),style='-';marker='.';else,style='--';marker='o';end
    if strcmp(view,'separation'),style='none';end
    h=line('Parent',ax,'XData',x,'YData',y,'LineStyle',style,'Marker',marker,'LineWidth',1.2);
    if strcmp(row.state,'accepted'),accepted(end+1,1)=h;else,unresolved(end+1,1)=h;end %#ok<AGROW>
end
handles=gobjects(0,1);labels={};
if ~isempty(accepted),handles(end+1)=accepted(1);labels{end+1}='stored accepted cells';end
if ~isempty(unresolved),handles(end+1)=unresolved(1);labels{end+1}='unresolved descriptors';end
if ~isempty(handles),legend(ax,handles,labels,'Interpreter','none','Location','best');end
if strcmp(view,'coverage')
    if strcmp(D.testVariable,'s'),xlabel(ax,'frequency omega');else,xlabel(ax,'forward unit-circle angle (rad)');end
    ylabel(ax,'stored record index');
else
    xlabel(ax,'stored record index');
    if isempty(units),ylabel(ax,'no stored finite metric');else,ylabel(ax,units{1},'Interpreter','none');end
end
grid(ax,'on');omitted=n-numel(idx);
heading=sprintf('Solver reported %s | %d retained / %g accepted cells',D.status,D.storedAcceptedCount,D.reportedAcceptedCount);
note=sprintf('%d unresolved descriptors; %d rows omitted; %d selected rows have no coordinate.', ...
    D.unresolvedDescriptorCount,omitted,noCoordinate);
if strcmp(view,'coverage'),meaning='Frequency projections do not independently prove parameter coverage.';
else,meaning='Stored separation/slack is not a spectral or parameter robustness margin.';end
if ~D.storageComplete,meaning=[meaning,' Accepted-leaf storage is incomplete or was not requested.'];end
metricCaption='';
if strcmp(view,'separation')
    for j=idx(:).'
        if isfinite(D.rows(j).metric),metricCaption=D.rows(j).metricName;break;end
    end
end
if isempty(metricCaption),title(ax,{heading,note,meaning},'Interpreter','none');
else,title(ax,{heading,note,metricCaption,meaning},'Interpreter','none');end
annotation=gobjects(0,1);
if isempty(handles)
    annotation=text(ax,.5,.5,'No plottable retained cells. Use StoreCertificates=true when computing R.', ...
        'Units','normalized','HorizontalAlignment','center','Interpreter','none');
end
H=struct('figure',fig,'axes',ax,'accepted',accepted,'unresolved',unresolved, ...
    'annotation',annotation,'view',view,'rowIndices',idx(:).','omittedRowCount',omitted, ...
    'selectedWithoutCoordinate',noCoordinate,'isCertificate',false,'data',D);
end
function tf=localCount(x)
tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==fix(x);
end
