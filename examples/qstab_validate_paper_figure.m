function qstab_validate_paper_figure(I)
%QSTAB_VALIDATE_PAPER_FIGURE Validate ordinary plotting data before side effects.
fields={'schema','id','article','sourceRevision','sourceAnchor','title','subtitle', ...
 'xlabel','ylabel','zlabel','dimension','equalAxes','isCertificate','curves'};
if ~isstruct(I)||~isscalar(I)||~all(isfield(I,fields))||~strcmp(I.schema,'qstab.paper-figure.1')
 error('qstab:paperFigure:Data','A scalar paper-figure data item is required.');
end
texts={'id','article','sourceRevision','sourceAnchor','title','subtitle','xlabel','ylabel','zlabel'};
for j=1:numel(texts)
 x=I.(texts{j});if ~ischar(x)||(~isrow(x)&&~isempty(x)),error('qstab:paperFigure:Text','Figure labels must be character vectors.');end
end
if isempty(regexp(I.id,'^[A-Za-z][A-Za-z0-9_]*$','once'))|| ...
 ~isequal(I.isCertificate,false)||~islogical(I.equalAxes)||~isscalar(I.equalAxes)|| ...
 ~isnumeric(I.dimension)||~isscalar(I.dimension)||~ismember(I.dimension,[2 3])
 error('qstab:paperFigure:Metadata','Invalid figure ID, dimensionality or assurance metadata.');
end
if ~isstruct(I.curves)||isempty(I.curves)||~isvector(I.curves)|| ...
 ~all(isfield(I.curves,{'x','y','z','label','role','lineStyle','marker','showLegend'}))
 error('qstab:paperFigure:Curves','Nonempty curve data required.');
end
roles={'computedSamples','analyticBoundary','outerBoundary','reference','guide','point','criticalPoint','control'};
for j=1:numel(I.curves)
 c=I.curves(j);
 for name={'x','y','z'}
  v=c.(name{1});if ~isnumeric(v)||~isreal(v)||(~isvector(v)&&~isempty(v))||any(~isfinite(v(:)))
   error('qstab:paperFigure:Coordinates','Coordinates must be finite real numeric vectors.');end
 end
 if isempty(c.x)||numel(c.x)~=numel(c.y)||(I.dimension==3&&numel(c.z)~=numel(c.x))||(I.dimension==2&&~isempty(c.z))
  error('qstab:paperFigure:Shape','Coordinate lengths must agree with figure dimensionality.');end
 if ~ischar(c.label)||(~isrow(c.label)&&~isempty(c.label))||~ischar(c.role)||~ismember(c.role,roles)|| ...
  ~ischar(c.lineStyle)||~ismember(c.lineStyle,{'-','--',':','-.','none'})|| ...
  ~ischar(c.marker)||~ismember(c.marker,{'none','o','x','+','s','d','.'})|| ...
  ~islogical(c.showLegend)||~isscalar(c.showLegend)
  error('qstab:paperFigure:Style','Invalid label, role, style or legend metadata.');end
end
end
