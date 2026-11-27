function files=plot_qstab_oneparameter(report,varargin)
%PLOT_QSTAB_ONEPARAMETER Plot computed scalar one-parameter article results.
%   FILES=plot_qstab_oneparameter(R,'Directory',folder) saves the E1 frequency
%   comparison and E3 nearest-gain geometry when those cases are present in R.
%   PNG, PDF, MATLAB FIG and accompanying CSV data are written. FILES gives
%   their absolute paths. 'Visible','on' also keeps the new figures open.
%   Existing figures and axes are not modified.
%
%   E1 uses the computed curves in R, including the restricted-slice quotient.
%   Black dots evaluate the reported reference frequencies; their placement
%   is not by itself a new certificate of global minimization.
%   E3 distinguishes the analytical sphere from the seven computed witness
%   checks. Samples alone do not prove the complete sphere characterization.
%   See also run_qstab_oneparameter.
ip=inputParser;ip.FunctionName='plot_qstab_oneparameter';ip.PartialMatching=false;
addParameter(ip,'Directory','',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'Visible','off',@(x)ischar(x)&&any(strcmp(x,{'on','off'})));
parse(ip,varargin{:});o=ip.Results;folder=char(o.Directory);
if isempty(folder),folder=report.outputDirectory;end
if ~isfolder(folder),mkdir(folder);end
files={};idx=find(strcmp({report.cases.id},'E1'),1);
if ~isempty(idx)&&isfield(report.cases(idx).output,'curves')
 D=report.cases(idx).output;G=D.curves;
 f=figure('Color','w','Visible',o.Visible,'Name','qstab E1: calibration margins','Position',[100 100 800 500]);
 c=onCleanup(@()closeIfHidden(f,o.Visible));
 ax=axes('Parent',f,'FontSize',11);hold(ax,'on');
 h(1)=plot(ax,G.omega,G.left,'Color',[.08 .34 .65],'LineWidth',1.7);
 h(2)=plot(ax,G.omega,G.right,'Color',[.76 .25 .15],'LineWidth',1.7);
 h(3)=plot(ax,G.omega,G.slice,'--','Color',[.25 .25 .25],'LineWidth',1.3);
 for k=1:numel(D.problems)
  Q=D.problems{k};scatter(ax,Q.frequency,Q.frequencyRadius,45,'k','filled','HandleVisibility','off');
 end
 grid(ax,'on');box(ax,'on');xlabel(ax,'Frequency \omega');ylabel(ax,'Minimum parameter modulus');
 title(ax,'E1: actuator and sensor calibration','FontSize',13);
 legend(ax,h,{'Left multiplication','Right multiplication','Fixed-slice quotient'},'Location','northwest');
 xlim(ax,[min(G.omega),max(G.omega)]);
 files=[files,saveFigure(f,folder,'E1_calibration_margins')]; %#ok<AGROW>
 p=fullfile(folder,'E1_calibration_margins.csv');writeCSV(p,'omega,left_radius,right_radius,fixed_slice_quotient',[G.omega(:),G.left(:),G.right(:),G.slice(:)]);files{end+1}=p;
 clear c
end
idx=find(strcmp({report.cases.id},'E3'),1);
if ~isempty(idx)&&isfield(report.cases(idx).output,'sphereParameters')
 D=report.cases(idx).output;T=D.sphereParameters;
 f=figure('Color','w','Visible',o.Visible,'Name','qstab E3: nearest quaternion gains','Position',[100 100 720 560]);
 c=onCleanup(@()closeIfHidden(f,o.Visible));ax=axes('Parent',f,'FontSize',11);hold(ax,'on');
 [x,y,z]=sphere(18);r=D.sphereRadius;
 hs=mesh(ax,r*x,r*y,r*z,'FaceColor','none','EdgeColor',[.55 .74 .82],'LineWidth',.4);
 hp=plot3(ax,T(2,:),T(3,:),T(4,:),'o','LineStyle','none','Color',[.72 .19 .12], ...
  'MarkerFaceColor',[.72 .19 .12],'MarkerSize',6);
 axis(ax,'equal');grid(ax,'on');box(ax,'on');view(ax,37,23);
 set(ax,'XTick',[-.5 0 .5],'YTick',[-.5 0 .5],'ZTick',[-.5 0 .5]);
 xlabel(ax,'i component');ylabel(ax,'j component');zlabel(ax,'k component');
 title(ax,'E3: nearest gains in the hyperplane Re(t) = -1/2','FontSize',13);
 set(ax,'Position',[.14 .22 .72 .65]);
 annotation(f,'textbox',[.08 .035 .84 .07],'String', ...
  'Blue mesh: analytical two-sphere. Red markers: computed witnesses.', ...
  'EdgeColor','none','HorizontalAlignment','center','FontSize',10);
 files=[files,saveFigure(f,folder,'E3_nearest_gain_sphere')]; %#ok<AGROW>
 p=fullfile(folder,'E3_nearest_gain_sphere.csv');writeCSV(p,'real,i,j,k,original_polynomial_residual',[T.',D.sphereResiduals(:)]);files{end+1}=p;
 clear c
end
end
function files=saveFigure(f,folder,name)
files={fullfile(folder,[name '.png']),fullfile(folder,[name '.pdf']),fullfile(folder,[name '.fig'])};
set(f,'PaperPositionMode','auto');
if exist('exportgraphics','file')==2
 exportgraphics(f,files{1},'Resolution',180);exportgraphics(f,files{2},'ContentType','vector');
else
 set(f,'PaperUnits','inches','PaperPosition',[0 0 8 5.5],'PaperSize',[8 5.5]);
 print(f,files{1},'-dpng','-r150');print(f,files{2},'-dpdf');
end
savefig(f,files{3});
end
function closeIfHidden(f,visible),if strcmp(visible,'off')&&ishghandle(f),close(f);end;end
function writeCSV(file,header,A)
f=fopen(file,'w');if f<0,error('qstab:oneparameterPlot:Write','Cannot write %s.',file);end
c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',header);fmt=[repmat('%.17g,',1,size(A,2)-1) '%.17g\n'];fprintf(f,fmt,A.');
end
