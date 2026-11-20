function paths=plot_qstab_intervals(report,varargin)
%PLOT_QSTAB_INTERVALS Export numerical illustrations from an existing report.
%   files = plot_qstab_intervals(R,'Directory',folder,'Visible','off')
%
%   R is a report returned by run_qstab_intervals. Only R.cases(...).output
%   is read: this function does not rerun stability tests, compute roots,
%   integrate trajectories, or execute a manual example. Missing, failed or
%   unselected cases are omitted. The complete report produces eight figures:
%     I03  running-square region and the three norm-value segments;
%     I06  Hurwitz square with a stable complete boundary;
%     I17  Schur square with a stable complete boundary;
%     I21  edge-interior reserve;
%     I19  five-gain gap and three time responses;
%     I20  the one-component frequency-margin curve.
%
%   Each figure is saved as PNG, vector PDF and editable MATLAB FIG. A CSV
%   file retains its plotted curve or geometric data. FILES is a cell column
%   of absolute paths, suitable for inclusion in the report ZIP. Directory
%   defaults to R.outputDirectory. Visible='off' closes only the new export
%   windows; 'on' leaves them available. Existing figures, axes and hold
%   states are preserved. No change is made to the input report or its data.
%
%   Region boundaries are drawn from the analytic formulas in the article
%   and supplement. Frequency curves and trajectories use computed report
%   arrays and are labelled as illustrations, not continuum certificates.
%
%   Usually request the figures through the main runner:
%     R = run_qstab_intervals('all','Plots',true);
%
%   See also run_qstab_intervals, qstab_interval_cases.
if ~(isstruct(report)&&isscalar(report)&&isfield(report,'cases')&&isstruct(report.cases))
    error('qstab:intervalGraphics:Report','Supply a report returned by run_qstab_intervals.');
end
if ~isempty(report.cases)&&(~isfield(report.cases,'id')||~isfield(report.cases,'output'))
    error('qstab:intervalGraphics:Report','Report cases must contain id and output fields.');
end
out='';if isfield(report,'outputDirectory'),out=report.outputDirectory;end
ip=inputParser;ip.FunctionName='plot_qstab_intervals';ip.PartialMatching=false;
addParameter(ip,'Directory',out,@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'Visible','off',@(x)ischar(x)&&any(strcmp(x,{'on','off'})));
parse(ip,varargin{:});o=ip.Results;out=char(o.Directory);
if isempty(out),error('qstab:intervalGraphics:Directory','Supply Directory or a report with outputDirectory.');end
if ~absolute(out),out=fullfile(pwd,out);end
if ~isfolder(out),[ok,msg]=mkdir(out);if ~ok,error('qstab:intervalGraphics:Directory','Cannot create output directory: %s',msg);end;end
oldFigure=get(0,'CurrentFigure');oldAxes=[];
if ~isempty(oldFigure)&&isgraphics(oldFigure),oldAxes=get(oldFigure,'CurrentAxes');end
restore=onCleanup(@()restoreCurrent(oldFigure,oldAxes)); %#ok<NASGU>
paths={};
blue=[.13 .36 .64];red=[.75 .19 .22];orange=[.85 .49 .12];gray=[.40 .43 .47];
lightBlue=[.88 .93 .98];lightRed=[.99 .90 .89];lightGray=[.94 .95 .96];

D=data('I03',{'rhos','radius','normValueRanges','edgeContacts'});
if ~isempty(D)
    rho=max(D.rhos);a=linspace(-rho,rho,501);b=sqrt(4+a.^2);
    [fig,ax]=canvas('I03_running_square');clean=onCleanup(@()finishFigure(fig,o.Visible)); %#ok<NASGU>
    fill(ax,[-rho rho rho -rho],[-rho -rho rho rho],lightRed,'EdgeColor','none');
    band=min(b,rho);fill(ax,[a fliplr(a)],[-band fliplr(band)],lightBlue,'EdgeColor','none');
    h1=plot(ax,a,b,'Color',gray,'LineWidth',1.4);plot(ax,a,-b,'Color',gray,'LineWidth',1.4);
    h2=plot(ax,[-rho rho rho -rho -rho],[-rho -rho rho rho -rho],'Color',blue,'LineWidth',1.5);
    h3=plot(ax,[-rho rho],[-rho rho],'--','Color',orange,'LineWidth',1.4);
    plot(ax,[-rho -rho rho rho],[-rho rho -rho rho],'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue,'MarkerSize',5);
    h4=plot(ax,0,rho,'^','MarkerFaceColor',red,'MarkerEdgeColor',red,'MarkerSize',7);
    h5=plot(ax,0,D.radius,'o','MarkerFaceColor','w','MarkerEdgeColor',gray,'MarkerSize',7,'LineWidth',1.4);
    plot(ax,D.edgeContacts,[rho rho],'o','MarkerFaceColor','w','MarkerEdgeColor',gray,'MarkerSize',5);
    text(ax,0,-.7,'Hurwitz','Color',blue,'HorizontalAlignment','center','FontSize',12);
    text(ax,0,2.65,'unstable','Color',red,'HorizontalAlignment','center','FontSize',10);
    axis(ax,'equal');xlim(ax,[-rho-.2 rho+.2]);ylim(ax,[-rho-.2 rho+.2]);
    xlabel(ax,'\alpha');ylabel(ax,'\beta');title(ax,'I03 | Running square: \beta^2 < 4 + \alpha^2');
    key(ax,[h1 h2 h3 h4 h5],{'Critical curve','Box','Correlated line','Unstable member','First contact'});
    note(fig,'Analytic regions; all four corners are stable. The full square and the correlated line are different families.');
    emit(fig,'I03_running_square',{'alpha','critical_beta_lower','critical_beta_upper'},[a(:),-b(:),b(:)]);clear clean

    [fig,ax]=canvas('I03_norm_segments');clean=onCleanup(@()finishFigure(fig,o.Visible)); %#ok<NASGU>
    yy=1:numel(D.rhos);xlo=min(D.normValueRanges(:,1));xhi=max(D.normValueRanges(:,2));
    plot(ax,[0 0],[.5 numel(yy)+.5],'--','Color',red,'LineWidth',1.3);
    for k=1:numel(yy)
        plot(ax,D.normValueRanges(k,:),[k k],'-o','Color',blue,'LineWidth',2.7, ...
            'MarkerFaceColor',blue,'MarkerSize',4);
    end
    [~,contactRow]=min(abs(D.rhos-D.radius));
    plot(ax,0,contactRow,'o','MarkerFaceColor','w','MarkerEdgeColor',red,'LineWidth',1.5,'MarkerSize',8);
    plot(ax,-4,numel(yy),'d','MarkerFaceColor',orange,'MarkerEdgeColor',orange,'MarkerSize',7);
    text(ax,-4,numel(yy)+.22,'vertices and correlated line','HorizontalAlignment','center','Color',orange,'FontSize',10);
    labels=cell(size(D.rhos));for k=1:numel(labels),labels{k}=sprintf('rho = %.3g',D.rhos(k));end
    set(ax,'YTick',yy,'YTickLabel',labels);ylim(ax,[.5 numel(yy)+.55]);xlim(ax,[xlo-1 xhi+1]);
    xlabel(ax,'Real part of N_p(i\omega)');
    % The rho tick labels and footer identify the rows without a long
    % redundant vertical label extending beyond some export backends.
    ylabel(ax,'');
    title(ax,'I03 | At \omega = 1: imaginary part is zero for every member');
    note(fig,'Exact whole-family segments. Vertical offsets identify different boxes; they are not imaginary coordinates.');
    emit(fig,'I03_norm_segments',{'rho','real_lower','real_upper','imaginary_value'}, ...
        [D.rhos(:),D.normValueRanges,zeros(numel(D.rhos),1)]);clear clean
end

D=data('I06',{'wholeValueRange','boundaryValueRange','boundaryMinimumDelta3'});
if ~isempty(D)
    th=linspace(0,2*pi,501);rr=sqrt(3);
    [fig,ax]=canvas('I06_hurwitz_boundary');clean=onCleanup(@()finishFigure(fig,o.Visible)); %#ok<NASGU>
    fill(ax,[-2.3 2.3 2.3 -2.3],[-2.3 -2.3 2.3 2.3],lightBlue,'EdgeColor','none');
    fill(ax,rr*cos(th),rr*sin(th),lightRed,'EdgeColor',gray,'LineWidth',1.3);
    h1=plot(ax,[-2 2 2 -2 -2],[-2 -2 2 2 -2],'Color',blue,'LineWidth',2);
    h2=plot(ax,0,0,'o','MarkerFaceColor',red,'MarkerEdgeColor',red,'MarkerSize',6);
    text(ax,0,.38,'y^2 + z^2 \leq 3','HorizontalAlignment','center','FontSize',12,'Color',red);
    text(ax,0,-.45,'not strictly Hurwitz','HorizontalAlignment','center','FontSize',10,'Color',red);
    axis(ax,'equal');xlim(ax,[-2.3 2.3]);ylim(ax,[-2.3 2.3]);xlabel(ax,'y');ylabel(ax,'z');
    title(ax,'I06 | A stable complete boundary around an unstable centre');
    key(ax,[h1 h2],{'All four continuous edges are stable','Unstable centre'});
    note(fig,sprintf('Analytic condition: y^2+z^2>3. At omega=1: full-square values [%g,%g], boundary [%g,%g].', ...
        D.wholeValueRange(1),D.wholeValueRange(2),D.boundaryValueRange(1),D.boundaryValueRange(2)));
    emit(fig,'I06_hurwitz_boundary',{'angle_rad','critical_y','critical_z'}, ...
        [th(:),rr*cos(th(:)),rr*sin(th(:))]);clear clean
end

D=data('I17',{'boundarySquaredNorm','completeBoundaryGaps'});
if ~isempty(D)
    th=linspace(0,2*pi,501);rin=sqrt(11/25);halfwidth=sqrt(min(D.boundarySquaredNorm));
    [fig,ax]=canvas('I17_schur_boundary');clean=onCleanup(@()finishFigure(fig,o.Visible)); %#ok<NASGU>
    fill(ax,[-1.15 1.15 1.15 -1.15],[-1.15 -1.15 1.15 1.15],lightGray,'EdgeColor','none');
    fill(ax,cos(th),sin(th),lightBlue,'EdgeColor',gray,'LineWidth',1.3);
    fill(ax,rin*cos(th),rin*sin(th),lightRed,'EdgeColor',gray,'LineWidth',1.3);
    h1=plot(ax,halfwidth*[-1 1 1 -1 -1],halfwidth*[-1 -1 1 1 -1],'Color',blue,'LineWidth',2);
    h2=plot(ax,0,0,'o','MarkerFaceColor',red,'MarkerEdgeColor',red,'MarkerSize',6);
    text(ax,0,.19,'t_1^2 + t_2^2 \leq 11/25','HorizontalAlignment','center','FontSize',11,'Color',red);
    text(ax,.85,0,'Schur','HorizontalAlignment','center','FontSize',10,'Color',blue);
    axis(ax,'equal');xlim(ax,[-1.15 1.15]);ylim(ax,[-1.15 1.15]);xlabel(ax,'t_1');ylabel(ax,'t_2');
    title(ax,'I17 | Monic Schur square: the stable region is an annulus');
    key(ax,[h1 h2],{'All four continuous edges are stable','Unstable centre'});
    note(fig,'Analytic condition: 11/25 < t_1^2+t_2^2 < 1. Both critical circles are excluded from strict stability.');
    emit(fig,'I17_schur_boundary',{'angle_rad','inner_t1','inner_t2','outer_t1','outer_t2'}, ...
        [th(:),rin*cos(th(:)),rin*sin(th(:)),cos(th(:)),sin(th(:))]);clear clean
end

D=data('I21',{'radius','rightRoot','vertexRayMinimum'});
if ~isempty(D)
    r=D.radius;th=linspace(0,2*pi,501);rr=sqrt(3);
    xx=[2-r 2+r 2+r 2-r 2-r];yy=2*r*[-1 -1 1 1 -1];
    [fig,ax]=canvas('I21_edge_reserve');clean=onCleanup(@()finishFigure(fig,o.Visible)); %#ok<NASGU>
    fill(ax,rr*cos(th),rr*sin(th),lightRed,'EdgeColor',gray,'LineWidth',1.2);
    h1=plot(ax,xx,yy,'Color',blue,'LineWidth',1.9);
    plot(ax,xx(1:4),yy(1:4),'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue,'MarkerSize',5);
    h2=plot(ax,2-r,0,'o','MarkerFaceColor',red,'MarkerEdgeColor',red,'MarkerSize',7);
    h3=plot(ax,2,0,'d','MarkerFaceColor','w','MarkerEdgeColor',orange,'LineWidth',1.5,'MarkerSize',6);
    text(ax,0,.3,'y^2+z^2 \leq 3','HorizontalAlignment','center','Color',red,'FontSize',11);
    axis(ax,'equal');xlim(ax,[-1.95 2.65]);ylim(ax,[-1.95 1.95]);xlabel(ax,'y');ylabel(ax,'z');
    title(ax,sprintf('I21 | Edge-interior first contact: r = %.10g',r));
    key(ax,[h1 h2 h3],{'Critical box','First contact','Nominal member'});
    note(fig,sprintf('All vertex rays satisfy y^2+z^2 >= %.6g > 3 for every r >= 0; the first contact is inside an edge.',D.vertexRayMinimum));
    emit(fig,'I21_edge_reserve',{'box_vertex_y','box_vertex_z'},[xx(:),yy(:)]);clear clean
end

D=data('I19',{'marginPolynomial','radius','budgets','time','trajectories','spectralAbscissae'});
if ~isempty(D)
    rho=linspace(0,1.25,401);gap=-polyval(D.marginPolynomial,rho)/1000;
    [fig,ax]=canvas('I19_regulation_gap');clean=onCleanup(@()finishFigure(fig,o.Visible)); %#ok<NASGU>
    h1=plot(ax,rho,gap,'Color',blue,'LineWidth',1.8);
    plot(ax,[0 1.25],[0 0],'-','Color',gray,'LineWidth',.8);
    h2=plot(ax,[1 1],[-.7 5.3],'--','Color',orange,'LineWidth',1.4);
    h3=plot(ax,D.radius,0,'o','MarkerFaceColor',red,'MarkerEdgeColor',red,'MarkerSize',7);
    xlim(ax,[0 1.25]);ylim(ax,[-.7 5.3]);xlabel(ax,'Tolerance multiplier \rho');ylabel(ax,'Worst-case Hurwitz gap F(\rho)');
    title(ax,sprintf('I19 | Five independent gains: reserve %.12g',D.radius));
    key(ax,[h1 h2 h3],{'Analytic whole-box gap','Design budget','First zero'});
    note(fig,'The gap formula holds on 0 <= rho <= 5/4. The displayed root is computed independently of this plotted grid.');
    emit(fig,'I19_regulation_gap',{'rho','worst_case_gap'},[rho(:),gap(:)]);clear clean

    [fig,ax]=canvas('I19_regulation_trajectories');clean=onCleanup(@()finishFigure(fig,o.Visible)); %#ok<NASGU>
    colors=[blue;orange;red];h=cell(1,numel(D.budgets));labels=cell(1,numel(D.budgets));
    for k=1:numel(D.budgets)
        h{k}=plot(ax,D.time,D.trajectories(:,k),'Color',colors(1+mod(k-1,3),:),'LineWidth',1.3); %#ok<AGROW>
        labels{k}=sprintf('rho = %.3g (abscissa %.5g /s)',D.budgets(k),D.spectralAbscissae(k));
    end
    set(ax,'YScale','log');xlim(ax,[min(D.time) max(D.time)]);ylim(ax,[1e-4 200]);
    xlabel(ax,'Time (s)');ylabel(ax,'Position-error norm (mm)');title(ax,'I19 | Illustrative responses from a 10 mm initial error');
    key(ax,[h{:}],labels);
    note(fig,'Computed trajectories at the gap-minimizing corners; these realizations need not maximize a transient norm.');
    headers=cell(1,1+numel(D.budgets));headers{1}='time_s';
    for k=1:numel(D.budgets),headers{k+1}=sprintf('error_norm_mm_rho_%g',D.budgets(k));end
    emit(fig,'I19_regulation_trajectories',headers,[D.time(:),D.trajectories]);clear clean
end

D=data('I20',{'frequency','frequencyMargin','radius','omega'});
if ~isempty(D)
    [fig,ax]=canvas('I20_two_gain_frequency');clean=onCleanup(@()finishFigure(fig,o.Visible)); %#ok<NASGU>
    h1=plot(ax,D.frequency,D.frequencyMargin,'Color',blue,'LineWidth',1.7);
    h2=plot(ax,[0 1.2],[D.radius D.radius],'--','Color',red,'LineWidth',1.2);
    plot(ax,D.omega,D.radius,'o','MarkerFaceColor',red,'MarkerEdgeColor',red,'MarkerSize',6);
    h3=plot(ax,[0 1.2],[1 1],':','Color',gray,'LineWidth',1.2);
    xlim(ax,[0 1.2]);ylim(ax,[0 18]);xlabel(ax,'Normalized frequency \omega');ylabel(ax,'Tolerance multiplier');
    title(ax,sprintf('I20 | One uncertain component: minimum %.12g',D.radius));
    key(ax,[h1 h2 h3],{'Two-point frequency formula','Algebraic reserve','Design budget'});
    note(fig,'The curve samples computed fixed-frequency margins. The global minimum follows from the separate cubic and contact witness.');
    emit(fig,'I20_two_gain_frequency',{'omega_normalized','frequency_margin'},[D.frequency(:),D.frequencyMargin(:)]);clear clean
end

    function D=data(id,fields)
        D=[];if isempty(report.cases),return;end
        j=find(strcmp({report.cases.id},id),1);
        if isempty(j),return;end
        R=report.cases(j);
        if isfield(R,'status')&&strcmp(R.status,'ERROR'),return;end
        if ~isstruct(R.output)||~isscalar(R.output)||~all(isfield(R.output,fields)),return;end
        D=R.output;
    end
    function [fig,ax]=canvas(name)
        fig=figure('Visible',o.Visible,'Name',name,'NumberTitle','off','Color','w', ...
            'Units','pixels','Position',[80 80 1000 680], ...
            'PaperUnits','inches','PaperPosition',[0 0 10 6.8],'PaperSize',[10 6.8]);
        ax=axes('Parent',fig,'Position',[.10 .28 .85 .61],'NextPlot','add', ...
            'FontSize',11,'Box','on','Layer','bottom','XGrid','on','YGrid','on','GridColor',[.88 .89 .91]);
    end
    function emit(fig,base,headers,values)
        names={fullfile(out,[base '.png']);fullfile(out,[base '.pdf']); ...
            fullfile(out,[base '.fig']);fullfile(out,[base '.csv'])};
        drawnow;
        print(fig,names{1},'-dpng','-r160');
        print(fig,names{2},'-dpdf','-painters');
        savefig(fig,names{3});
        writeCSV(names{4},headers,values);
        for j=1:numel(names)
            if ~isfile(names{j}),error('qstab:intervalGraphics:Export','Missing exported file %s.',names{j});end
        end
        paths=[paths;names];
    end
end
function key(ax,h,labels)
% Draw a compact key in fixed figure coordinates. This also avoids export
% backends repositioning a legend over the data or outside the page.
fig=ancestor(ax,'figure');n=numel(labels);columns=min(3,n);width=.88/columns;
for k=1:n
    row=floor((k-1)/columns);col=mod(k-1,columns);
    x=.08+col*width;y=.15-.035*row;
    a=axes('Parent',fig,'Position',[x y-.009 .035 .020], ...
        'Color','none','Visible','off','XLim',[0 1],'YLim',[0 1]);
    style=get(h(k),'LineStyle');marker=get(h(k),'Marker');
    if strcmp(style,'none'),xx=.5;else,xx=[.05 .95];end
    plot(a,xx,.5*ones(size(xx)),'Color',get(h(k),'Color'), ...
        'LineStyle',style,'LineWidth',max(1,get(h(k),'LineWidth')), ...
        'Marker',marker,'MarkerSize',get(h(k),'MarkerSize'), ...
        'MarkerFaceColor',get(h(k),'MarkerFaceColor'), ...
        'MarkerEdgeColor',get(h(k),'MarkerEdgeColor'));
    axis(a,'off');
    annotation(fig,'textbox',[x+.043 y-.016 width-.045 .035], ...
        'String',labels{k},'Interpreter','none','EdgeColor','none','FontSize',9, ...
        'VerticalAlignment','middle');
end
set(fig,'CurrentAxes',ax);
end
function note(fig,message)
words=strsplit(message);lines={''};
for k=1:numel(words)
    if numel(lines{end})+1+numel(words{k})>110,lines{end+1}=words{k}; %#ok<AGROW>
    elseif isempty(lines{end}),lines{end}=words{k};
    else,lines{end}=[lines{end} ' ' words{k}];end
end
previous=get(fig,'CurrentAxes');
a=axes('Parent',fig,'Position',[.06 .01 .91 .07],'Color','none', ...
    'Visible','off','XLim',[0 1],'YLim',[0 1]);
for k=1:numel(lines)
    text(a,.5,1-(k-.5)/numel(lines),lines{k},'Interpreter','none','FontSize',9, ...
        'HorizontalAlignment','center','VerticalAlignment','middle','Color',[.35 .37 .40]);
end
if ~isempty(previous)&&isgraphics(previous),set(fig,'CurrentAxes',previous);end
end
function finishFigure(fig,visible)
if strcmp(visible,'off')&&isgraphics(fig),close(fig);end
end
function restoreCurrent(fig,ax)
if ~isempty(fig)&&isgraphics(fig)
    set(0,'CurrentFigure',fig);
    if ~isempty(ax)&&isgraphics(ax),set(fig,'CurrentAxes',ax);end
elseif isempty(fig)
    set(0,'CurrentFigure',[]);
end
end
function y=absolute(x)
y=~isempty(x)&&(x(1)=='/'||x(1)=='\'||(numel(x)>1&&x(2)==':'));
end
function writeCSV(file,headers,A)
f=fopen(file,'w');if f<0,error('qstab:intervalGraphics:Write','Cannot write %s.',file);end
cleanup=onCleanup(@()fclose(f)); %#ok<NASGU>
for k=1:numel(headers),if k>1,fprintf(f,',');end;fprintf(f,'%s',headers{k});end
fprintf(f,'\n');
for j=1:size(A,1)
    for k=1:size(A,2),if k>1,fprintf(f,',');end;fprintf(f,'%.17g',A(j,k));end
    fprintf(f,'\n');
end
end
