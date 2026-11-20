function tests=qstab_cases_intervalgraphics()
%QSTAB_CASES_INTERVALGRAPHICS Report-only export and caller graphics state.
tests={ ...
 'All eight interval illustrations export PNG PDF FIG and unaltered numerical CSV data',@allKinds; ...
 'Selected report export preserves an existing axes, hold state and current figure',@preserveCaller; ...
 'Absent or failed cases do not invent plots and malformed reports fail before drawing',@absentCases ...
};
end

function allKinds()
R=fixture();copy=R;out=tempname;mkdir(out);cleanup=onCleanup(@()removeTree(out)); %#ok<NASGU>
old=findall(0,'Type','figure');files=plot_qstab_intervals(R,'Directory',out,'Visible','off');
assert(numel(files)==32&&all(cellfun(@isfile,files)));
for ext={'.png','.pdf','.fig','.csv'}
    names=dir(fullfile(out,['*' ext{1}]));assert(numel(names)==8&&all([names.bytes]>0));
end
D=dlmread(fullfile(out,'I19_regulation_trajectories.csv'),',',1,0);
I=R.cases(strcmp({R.cases.id},'I19')).output;
assert(isequal(D,[I.time,I.trajectories]));
D=dlmread(fullfile(out,'I20_two_gain_frequency.csv'),',',1,0);
I=R.cases(strcmp({R.cases.id},'I20')).output;
assert(isequal(D,[I.frequency(:),I.frequencyMargin(:)]));
assert(isequaln(R,copy)&&isequal(findall(0,'Type','figure'),old));
end

function preserveCaller()
out=tempname;mkdir(out);cleanup=onCleanup(@()removeTree(out)); %#ok<NASGU>
f=figure('Visible','off');cf=onCleanup(@()closeIfOpen(f)); %#ok<NASGU>
ax=axes('Parent',f);line=plot(ax,[0 1],[2 3]);set(ax,'NextPlot','replace');
old=findall(0,'Type','figure');oldpath=path;folder=pwd;state=rng;
R=fixture();R.cases=R.cases(strcmp({R.cases.id},'I20'));copy=R;
files=plot_qstab_intervals(R,'Directory',out,'Visible','off');
assert(numel(files)==4&&all(cellfun(@isfile,files)));
assert(isgraphics(line)&&isequal(get(line,'YData'),[2 3]));
assert(strcmp(get(ax,'NextPlot'),'replace')&&isequal(get(0,'CurrentFigure'),f));
assert(isequal(get(f,'CurrentAxes'),ax)&&isequal(findall(0,'Type','figure'),old));
assert(strcmp(path,oldpath)&&strcmp(pwd,folder)&&isequal(rng,state)&&isequaln(R,copy));
end

function absentCases()
out=tempname;mkdir(out);cleanup=onCleanup(@()removeTree(out)); %#ok<NASGU>
R=struct('cases',struct('id','I01','status','MATCHED','output',struct('norm',[2 2 1])));
old=findall(0,'Type','figure');files=plot_qstab_intervals(R,'Directory',out);
assert(isempty(files)&&isequal(findall(0,'Type','figure'),old));
R=fixture();R.cases=R.cases(1);R.cases.status='ERROR';
assert(isempty(plot_qstab_intervals(R,'Directory',out)));
R=struct('cases',struct([]));assert(isempty(plot_qstab_intervals(R,'Directory',out)));
caught=false;try,plot_qstab_intervals(struct('cases',{{}}),'Directory',out);
catch err,caught=strcmp(err.identifier,'qstab:intervalGraphics:Report');end
assert(caught&&isequal(findall(0,'Type','figure'),old));
end

function R=fixture()
% Precomputed plotting fixtures: no numerical solver is invoked by this test.
R=struct('cases',struct('id',{},'status',{},'output',{}));
D=struct('rhos',[1.5 2 3],'radius',2,'normValueRanges',[-6.25 -1.75;-8 0;-13 5], ...
 'edgeContacts',[-sqrt(5) sqrt(5)]);put('I03',D);
D=struct('wholeValueRange',[-5 3],'boundaryValueRange',[-5 -1],'boundaryMinimumDelta3',4);put('I06',D);
D=struct('boundarySquaredNorm',[.49 .98],'completeBoundaryGaps',[.05 .02]);put('I17',D);
D=struct('radius',2-sqrt(3),'rightRoot',[0;-.5;0;-sqrt(3)/2],'vertexRayMinimum',16/5);put('I21',D);
D=struct('marginPolynomial',[51 -920 5500 -5000],'radius',1.0987217684952284, ...
 'budgets',[0 1 1.2],'time',[0;1;2], ...
 'trajectories',[10 10 10;7.472291686 9.4 10.5;5.746815095 8.9 11.3], ...
 'spectralAbscissae',[-.2620761414007129 -.0524051008658224 .0658486838669244]);put('I19',D);
D=struct('frequency',[.15 .5 .6969070706851537 sqrt(.5) .9 1.2], ...
 'frequencyMargin',[120.298801905 25 5.143205348290382 10 10 28], ...
 'radius',5.143205348290382,'omega',.6969070706851537);put('I20',D);
    function put(id,D)
        R.cases(end+1)=struct('id',id,'status','MATCHED','output',D);
    end
end
function closeIfOpen(f)
if isgraphics(f),close(f);end
end
function removeTree(p)
if isfolder(p),rmdir(p,'s');end
end
