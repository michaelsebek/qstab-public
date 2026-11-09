function X=export_qstab_paper_figures(F,varargin)
%EXPORT_QSTAB_PAPER_FIGURES Export SAVED data without recomputing solvers.
% X=export_qstab_paper_figures(F) writes PDF/PNG/FIG/CSV and a text manifest
% into <package root>/qstab_paper_figures, independently of the current folder.
% Options: OutputDir, Formats={'pdf','png','fig','csv'}, Resolution=180,
% Overwrite=true. Relative OutputDir is resolved against this package root.
% Only the selected figure files and figures_manifest.txt are replaced.
% All figures are generated in a staging folder before committing any target.
% No test reports/ZIPs are read, written or renamed. No figure handles are
% stored in the returned metadata. CSV z-coordinates are zero for 2-D plots.
root=fileparts(fileparts(mfilename('fullpath')));
ip=inputParser;ip.FunctionName='export_qstab_paper_figures';ip.PartialMatching=false;
addParameter(ip,'OutputDir',fullfile(root,'qstab_paper_figures'));
addParameter(ip,'Formats',{'pdf','png','fig','csv'});
addParameter(ip,'Resolution',180,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x<=2400&&x==fix(x));
addParameter(ip,'Overwrite',true,@(x)islogical(x)&&isscalar(x));parse(ip,varargin{:});o=ip.Results;
if ~isstruct(F)||~isscalar(F)||~all(isfield(F,{'schema','items','ok','toolbox','runtime','nativeQuaternionInputs'}))|| ...
 ~strcmp(F.schema,'qstab.paper-figures.1')||~isequal(F.ok,true)||~isstruct(F.items)||isempty(F.items)
 error('qstab:paperExport:Data','A checked paper-figure data package is required.');end
for j=1:numel(F.items),qstab_validate_paper_figure(F.items(j));end
if numel(unique({F.items.id}))~=numel(F.items),error('qstab:paperExport:Data','Duplicate figure IDs.');end
formats=o.Formats;
if ischar(formats)&&isrow(formats)&&~isempty(formats),formats={formats};elseif isstring(formats)&&isvector(formats),formats=cellstr(formats);end
if ~iscell(formats)||~isvector(formats)||isempty(formats),error('qstab:paperExport:Formats','Formats must be a nonempty cell vector.');end
for j=1:numel(formats)
 if isstring(formats{j})&&isscalar(formats{j}),formats{j}=char(formats{j});end
 if ~ischar(formats{j})||~ismember(formats{j},{'pdf','png','fig','csv'}),error('qstab:paperExport:Formats','Supported formats are pdf, png, fig, csv.');end
end
if numel(unique(formats))~=numel(formats),error('qstab:paperExport:Formats','Duplicate formats.');end
out=o.OutputDir;if isstring(out)&&isscalar(out),out=char(out);end
if ~ischar(out)||~isrow(out)||isempty(strtrim(out))||any(out==char(0))
 error('qstab:paperExport:OutputDir','OutputDir must be a nonempty path.');end
if ispc
 absolute=~isempty(regexp(out,'^[A-Za-z]:[\\/]','once'))|| ...
  (numel(out)>=2&&((out(1)=='\'&&out(2)=='\')||(out(1)=='/'&&out(2)=='/')));
 if ~absolute&&(~isempty(regexp(out,'^[A-Za-z]:','once'))||out(1)=='\'||out(1)=='/')
  error('qstab:paperExport:OutputDir','Use an absolute drive/UNC path or a package-relative path.');end
else,absolute=out(1)=='/';end
if ~absolute,out=fullfile(root,out);end
if exist(out,'file')==2,error('qstab:paperExport:OutputDir','OutputDir is a file.');end
needGraphics=any(~strcmp(formats,'csv'));
if any(ismember(formats,{'png','pdf'}))&&exist('exportgraphics','file')==0&&exist('exportgraphics','builtin')==0
 error('qstab:paperExport:Unavailable','PNG/PDF export requires MATLAB exportgraphics (R2020a or later).');end
names={'figures_manifest.txt'};
for j=1:numel(F.items),for k=1:numel(formats),names{end+1}=[F.items(j).id,'.',formats{k}];end,end
if ~o.Overwrite
 for j=1:numel(names),if exist(fullfile(out,names{j}),'file')~=0,error('qstab:paperExport:Exists','Refusing to overwrite %s.',names{j});end,end
end
if exist(out,'dir')~=7,[ok,msg]=mkdir(out);if ~ok,error('qstab:paperExport:Directory','%s',msg);end,end
stage=tempname(out);[ok,msg]=mkdir(stage);if ~ok,error('qstab:paperExport:Directory','%s',msg);end
cleanStage=onCleanup(@()removeStage(stage)); %#ok<NASGU>
oldFig=get(0,'CurrentFigure');restoreCurrent=onCleanup(@()restoreFigure(oldFig)); %#ok<NASGU>
X=struct('schema','qstab.paper-export.1','status','RUNNING','ok',false, ...
 'outputDirectory',out,'formats',{formats},'files',struct('id',{},'format',{},'file',{},'bytes',{}), ...
 'manifestFile',fullfile(out,'figures_manifest.txt'),'sourceToolbox',F.toolbox, ...
 'dataRuntime',F.runtime,'nativeQuaternionData',F.nativeQuaternionInputs, ...
 'exportRuntime',version,'createdAt',datestr(now,30),'solverCalls',0);
for j=1:numel(F.items)
 exportOne(F.items(j),stage,formats,o.Resolution,needGraphics);
 for k=1:numel(formats)
  n=[F.items(j).id,'.',formats{k}];p=fullfile(stage,n);v=dir(p);
  if isempty(v)||v.bytes==0,error('qstab:paperExport:Empty','Empty export %s.',n);end
  checkSignature(p,formats{k});
  X.files(end+1)=struct('id',F.items(j).id,'format',formats{k},'file',fullfile(out,n),'bytes',v.bytes);
 end
end
writeManifest(fullfile(stage,'figures_manifest.txt'),F,X);
for j=1:numel(names)
 [ok,msg]=movefile(fullfile(stage,names{j}),fullfile(out,names{j}),'f');
 if ~ok,error('qstab:paperExport:Commit','Export commit failed: %s. Some earlier files may have been committed.',msg);end
end
X.status='COMPLETED';X.ok=true;
end
function exportOne(I,stage,formats,resolution,needGraphics)
H=[];
if needGraphics,H=plot_qstab_paper_figure(I,'Visible','off');end
cleanup=onCleanup(@()closePlot(H)); %#ok<NASGU>
for j=1:numel(formats)
 file=fullfile(stage,[I.id,'.',formats{j}]);
 switch formats{j}
  case 'pdf',drawnow;exportgraphics(H.figure,file,'ContentType','vector');
  case 'png',drawnow;exportgraphics(H.figure,file,'Resolution',resolution);
  case 'fig',savefig(H.figure,file);
  case 'csv',writeCSV(file,I);
 end
end
end
function writeCSV(file,I)
fid=fopen(file,'w');if fid<0,error('qstab:paperExport:Write','Cannot write %s.',file);end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'curve,point,x,y,z\n');
for j=1:numel(I.curves)
 c=I.curves(j);z=c.z;if isempty(z),z=zeros(size(c.x));end
 for k=1:numel(c.x),fprintf(fid,'%d,%d,%.17g,%.17g,%.17g\n',j,k,c.x(k),c.y(k),z(k));end
end
end
function writeManifest(file,F,X)
fid=fopen(file,'w');if fid<0,error('qstab:paperExport:Write','Cannot write manifest.');end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'qstab selected mathematical illustrations\nSource data version: %s\nData runtime: %s\nExport runtime: %s\n',F.toolbox.version,F.runtime,X.exportRuntime);
fprintf(fid,'No numerical solver was called during export. No full-article reproduction or interval-arithmetic certification is asserted.\n');
fprintf(fid,'CSV columns: curve,point,x,y,z. z=0 for 2-D data. Coordinates are written with 17 significant digits.\n');
for j=1:numel(F.items)
 I=F.items(j);fprintf(fid,'\n%s | %s | %s | %s\n%s\n%s\n',I.id,I.article,I.sourceRevision,I.sourceAnchor,I.title,I.subtitle);
 for k=1:numel(I.curves),fprintf(fid,'  curve %d | %s | %s\n',k,I.curves(k).role,I.curves(k).label);end
end
for j=1:numel(X.files),f=X.files(j);fprintf(fid,'\n%s.%s: %d bytes',f.id,f.format,f.bytes);end
fprintf(fid,'\n');
end
function checkSignature(file,fmt)
fid=fopen(file,'rb');if fid<0,error('qstab:paperExport:Read','Cannot read generated file.');end
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
h=fread(fid,8,'*uint8').';
if strcmp(fmt,'pdf')&&(numel(h)<5||~isequal(char(h(1:5)),'%PDF-')),error('qstab:paperExport:Signature','Invalid PDF signature.');end
if strcmp(fmt,'png')&&~isequal(h,uint8([137 80 78 71 13 10 26 10])),error('qstab:paperExport:Signature','Invalid PNG signature.');end
end
function closePlot(H),if isstruct(H)&&isfield(H,'figure')&&isgraphics(H.figure),close(H.figure);end,end
function removeStage(p),if exist(p,'dir')==7,rmdir(p,'s');end,end
function restoreFigure(f)
if isempty(f),set(0,'CurrentFigure',[]);elseif isgraphics(f,'figure'),set(0,'CurrentFigure',f);end
end
