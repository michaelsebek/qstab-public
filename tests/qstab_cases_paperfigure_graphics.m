function T=qstab_cases_paperfigure_graphics()
% Hidden graphics; actual export smoke test isolated in a temporary folder.
T={ ...
 'All eight saved paper figures render without changing their data',@allFigures; ...
 'Saved data are the exact line coordinates used by the paper renderer',@coordinates; ...
 'Existing axes children and NextPlot state survive paper-figure rendering',@parent; ...
 'Malformed saved coordinates are rejected before a new figure is created',@invalid; ...
 'Bernstein local hulls and unattainable controls retain different graphics roles',@roles; ...
 'Tetrahedron uses three actual parameter coordinates, not a spectral stability projection',@threeDimensions; ...
 'Paper PDF PNG FIG and CSV export creates readable outputs and closes its own figures',@exportSmoke};
end
function F=go(id),F=qstab_paper_figures('Figures',id,'UseQuaternion',false,'Samples',21);assert(F.ok);end
function allFigures()
F=go({});for j=1:numel(F.items),renderOne(F.items(j));end
end
function renderOne(I)
B=I;H=plot_qstab_paper_figure(I,'Visible','off');c=onCleanup(@()close(H.figure)); %#ok<NASGU>
assert(numel(H.lines)==numel(I.curves)&&isequaln(I,B)&&~H.isCertificate);
end
function coordinates()
F=go('qb01_sweep');H=plot_qstab_paper_figure(F.items,'Visible','off');c=onCleanup(@()close(H.figure)); %#ok<NASGU>
for j=1:numel(H.lines)
 x=get(H.lines{j},'XData');y=get(H.lines{j},'YData');
 assert(isequal(x(:),F.items.curves(j).x)&&isequal(y(:),F.items.curves(j).y));
end
end
function parent()
f=figure('Visible','off');c=onCleanup(@()close(f));ax=axes('Parent',f); %#ok<NASGU>
old=line('Parent',ax,'XData',[0 1],'YData',[0 1]);set(ax,'NextPlot','replacechildren');F=go('qb02_hull');
H=plot_qstab_paper_figure(F.items,'Parent',ax);assert(isgraphics(old)&&H.figure==f&&strcmp(get(ax,'NextPlot'),'replacechildren'));
end
function invalid()
F=go('qb02_counter');I=F.items;I.curves(1).y(1)=Inf;n=numel(findall(0,'Type','figure'));caught=false;
try,plot_qstab_paper_figure(I,'Visible','off');catch,caught=true;end
assert(caught&&numel(findall(0,'Type','figure'))==n);
end
function roles()
F=go('qb03_parabola');H=plot_qstab_paper_figure(F.items,'Visible','off');c=onCleanup(@()close(H.figure)); %#ok<NASGU>
assert(numel(findobj(H.axes,'Tag','qstab.paper.outerBoundary'))==3&&numel(findobj(H.axes,'Tag','qstab.paper.control'))==1);
end
function threeDimensions()
F=go('qb03_tetrahedron');H=plot_qstab_paper_figure(F.items,'Visible','off');c=onCleanup(@()close(H.figure)); %#ok<NASGU>
h=findobj(H.axes,'Tag','qstab.paper.criticalPoint');assert(numel(h)==1&&get(h,'XData')==0&&get(h,'YData')==0&&get(h,'ZData')==0);
assert(strcmp(get(get(H.axes,'ZLabel'),'String'),'y_k'));
end
function exportSmoke()
F=go('qb03_parabola');p=tempname;mkdir(p);clean=onCleanup(@()rmdir(p,'s')); %#ok<NASGU>
old=findall(0,'Type','figure');X=export_qstab_paper_figures(F,'OutputDir',p);
assert(X.ok&&numel(X.files)==4&&isequal(old,findall(0,'Type','figure'))&&X.solverCalls==0);
I=imread(fullfile(p,'qb03_parabola.png'));assert(~isempty(I)&&size(I,1)>50&&size(I,2)>50);
f=openfig(fullfile(p,'qb03_parabola.fig'),'new','invisible');c=onCleanup(@()close(f)); %#ok<NASGU>
assert(numel(findobj(f,'Tag','qstab.paper.outerBoundary'))==3);
for j=1:numel(X.files),assert(exist(X.files(j).file,'file')==2&&X.files(j).bytes>0);end
end
