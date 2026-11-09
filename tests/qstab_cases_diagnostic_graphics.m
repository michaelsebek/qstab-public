function tests = qstab_cases_diagnostic_graphics()
%QSTAB_CASES_DIAGNOSTIC_GRAPHICS Hidden-figure checks, no exported files.
tests={ ...
 'Outer polygons, controls, admissible samples and zero have separate objects',@enclosureObjects; ...
 'Local enclosure union is not replaced by one global hull',@unionObjects; ...
 'Enclosure plots preserve existing objects and exact axes NextPlot state',@existingAxes; ...
 'Mismatched sample arguments fail before creating a figure',@badSamples; ...
 'Coverage plots retain unresolved obligations and explicit display truncation',@coverage; ...
 'Lyapunov separation plots do not relabel slack as a stability margin',@slack; ...
 'A stable report without stored cells has an explicit empty-data annotation',@emptyStorage; ...
 'Delay witness plots keep input and forward coordinates distinct',@delayWitness; ...
 'Degree loss produces a diagnosis without invented finite roots',@degreeLoss};
end
function enclosureObjects()
V={[1 1],[3 1]};[~,E] = qstab.valueenclosure(V,1i);[~,S] = qstab.valueset(V,1i,'Weights',[.5 .5]);
[~,H] = qstab.plotEnclosure(E,'Samples',S,'Visible','off');c=onCleanup(@()close(H.figure)); %#ok<NASGU>
assert(numel(H.outlines)==1&&numel(H.controls)==1&&~H.isCertificate);
assert(isgraphics(H.samples,'line')&&abs(get(H.samples,'XData')-3)<1e-12);
assert(H.displaysAllRegions&&isgraphics(H.zero,'line'));
end
function unionObjects()
V={P(-2),P(2)};[~,E] = qstab.valueenclosure(V,0,'VertexCells',{[1 .5;0 .5],[.5 0;.5 1]});
[~,H] = qstab.plotEnclosure(E,'Visible','off');c=onCleanup(@()close(H.figure)); %#ok<NASGU>
assert(numel(H.outlines)==2&&H.displaysAllRegions);
for j=1:2,assert(min(get(H.outlines(j),'XData'))>=1-1e-12);end
end
function existingAxes()
f=figure('Visible','off');c=onCleanup(@()close(f)); %#ok<NASGU>
ax=axes('Parent',f);old=line('Parent',ax,'XData',[0 1],'YData',[0 1]);set(ax,'NextPlot','replacechildren');
[~,E] = qstab.valueenclosure({[1 1],[2 1]},1i);qstab.plotEnclosure(E,'Parent',ax);
assert(isgraphics(old)&&strcmp(get(ax,'NextPlot'),'replacechildren'));
[~,R] = qstab.zeroexclude({[1 1],[2 1]},'StoreCertificates',true,'MaxTime',Inf);
qstab.plotCertificate(R,'Parent',ax);assert(isgraphics(old)&&strcmp(get(ax,'NextPlot'),'replacechildren'));
end
function badSamples()
[~,E] = qstab.valueenclosure({[1 1],[2 1]},1i);[~,S] = qstab.valueset([1 1],2i);
old=numel(findall(0,'Type','figure'));caught=false;
try,qstab.plotEnclosure(E,'Samples',S,'Visible','off');catch ME,caught=strcmp(ME.identifier,'qstab:plotEnclosure:Lambda');end
assert(caught&&numel(findall(0,'Type','figure'))==old);
end
function coverage()
V=repmat({[1 1]},1,5);[~,R] = qstab.zeroexclude(V,'Reduction','four-generator', ...
    'MaxSubfamilies',1,'StoreCertificates',true,'MaxTime',Inf);
[~,D] = qstab.certificateData(R);[~,H] = qstab.plotCertificate(D,'Visible','off','MaxRows',1);
c=onCleanup(@()close(H.figure)); %#ok<NASGU>
assert(H.omittedRowCount==numel(D.rows)-1&&~H.isCertificate&&strcmp(H.data.status,'undetermined'));
end
function slack()
[~,R] = qstab.robuststable({[1 1]},'StoreCertificates',true,'MaxTime',Inf);
[~,H] = qstab.plotCertificate(R,'Visible','off');c=onCleanup(@()close(H.figure)); %#ok<NASGU>
assert(strcmp(H.view,'separation'));
assert(strcmp(get(get(H.axes,'YLabel'),'String'),'normalizedLyapunovInequality'));
end
function emptyStorage()
[~,R] = qstab.zeroexclude({[1 1],[2 1]},'MaxTime',Inf);
[~,H] = qstab.plotCertificate(R,'Visible','off');c=onCleanup(@()close(H.figure)); %#ok<NASGU>
assert(isgraphics(H.annotation,'text')&&strcmp(H.data.status,'stable')&&~H.data.storageComplete);
end
function delayWitness()
V={[1 -1.2],[1 -1.4]};[~,R] = qstab.zeroexclude(V,'Variable','zinv','MaxTime',Inf);[~,W] = qstab.witnessData(V,R);
[~,H] = qstab.plotWitness(W,'Visible','off');c=onCleanup(@()close(H.figure)); %#ok<NASGU>
assert(strcmp(H.variable,'zinv')&&abs(H.canonicalLambda-1/1.3)<1e-10);
[~,G] = qstab.plotWitness(W,'Representation','forward','Dimension',3,'Visible','off');d=onCleanup(@()close(G.figure)); %#ok<NASGU>
assert(strcmp(G.variable,'z')&&abs(G.canonicalLambda-1.3)<1e-10);
end
function degreeLoss()
V={[1 1],[-1 -1]};[~,R] = qstab.zeroexclude(V,'MaxTime',Inf);[~,W] = qstab.witnessData(V,R);
[~,H] = qstab.plotWitness(W,'Visible','off');c=onCleanup(@()close(H.figure)); %#ok<NASGU>
assert(isgraphics(H.annotation,'text')&&W.isDegreeLoss&&isempty(W.inputRepresentatives));
W.isAdmissible=false;old=numel(findall(0,'Type','figure'));caught=false;
try,qstab.plotWitness(W,'Visible','off');catch ME,caught=strcmp(ME.identifier,'qstab:plotWitness:Witness');end
assert(caught&&numel(findall(0,'Type','figure'))==old);
end
function x=P(t)
x=struct('a',[1 1],'b',[0 0],'c',[t 0],'d',[0 0]);
end
