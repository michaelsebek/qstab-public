function [out,figures] = run_qstab_diagnostics(varargin)
%RUN_QSTAB_DIAGNOSTICS Reproducible enclosure, coverage and witness examples.
% OUT=run_qstab_diagnostics('Plots',false) returns ONLY ordinary data.
% [OUT,FIGURES]=run_qstab_diagnostics('Plots',true) also draws eight separate
% figures. Handles are returned ONLY as the second output, never inside OUT.
% No path, directory, RNG or archive modifications are made. test_qstab can
% retain OUT in report.diagnosticExamples in the existing two-entry ZIP.
% This is a small diagnostic study, not a rerun of the timing benchmark.
ip=inputParser;addParameter(ip,'Plots',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Visible','on',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'Verbose',true,@(x)islogical(x)&&isscalar(x));parse(ip,varargin{:});o=ip.Results;
figures=[];if o.Plots,figures=gobjects(0,1);end
out=struct('schema','qstab.diagnostic-examples.1','ok',false,'sourceVersion',qstab.about(), ...
    'assurance','ordinaryDoubleDiagnostics; no new validated certificates');
V={poly([1 1;0 0;-2 0;0 0]),poly([1 1;0 0;2 0;0 0])};
t=linspace(0,1,41).';cells={[1 .5;0 .5],[.5 0;.5 1]};
out.bernstein=struct('family',{V});
[~,out.bernstein.coarse]=qstab.valueenclosure(V,0);
[~,out.bernstein.refined]=qstab.valueenclosure(V,0,'VertexCells',cells);
[~,out.bernstein.samples]=qstab.valueset(V,0,'Weights',[1-t,t]);
[~,out.bernstein.result] = qstab.zeroexclude(V,'StoreCertificates',true,'MaxTime',Inf);
[~,out.bernstein.coverage] = qstab.certificateData(out.bernstein.result);
[~,out.bernstein.leaf] = qstab.certificateEnclosure(V,out.bernstein.result);
F=struct('center',poly([0 0;1 1;0 0;0 0]),'generators',{{[1 1]}},'bounds',[-2 2]);
out.correlated=struct('family',F);
[~,out.correlated.enclosure]=qstab.valueenclosure(F,1i);
[~,out.correlated.result] = qstab.zeroexclude(F,'StoreCertificates',true,'MaxTime',Inf);
[~,out.correlated.coverage] = qstab.certificateData(out.correlated.result);
[~,out.correlated.leaf] = qstab.certificateEnclosure([],out.correlated.result);
L=poly([1 1;-.05 0;-.05 0;-.05 0]);U=poly([1.2 1;.05 0;.05 0;.05 0]);
out.interval=struct('family',{{L,U}});
[~,out.interval.result] = qstab.zeroexclude({L,U},'Uncertainty','interval','StoreCertificates',true,'MaxTime',Inf);
[~,out.interval.coverage] = qstab.certificateData(out.interval.result);
[~,out.interval.leaf] = qstab.certificateEnclosure([],out.interval.result);
P=repmat({[1 1]},1,5);out.budget=struct('family',{P});
[~,out.budget.result] = qstab.zeroexclude(P,'Reduction','four-generator','MaxSubfamilies',1, ...
    'StoreCertificates',true,'MaxTime',Inf);[~,out.budget.coverage] = qstab.certificateData(out.budget.result);
Y=2*[1 1 1;1 -1 -1;-1 1 -1;-1 -1 1];T=cell(1,4);
for j=1:4,T{j}=poly([-1 1 1;Y(j,1) 3*Y(j,1) 0;Y(j,2) 3*Y(j,2) 0;Y(j,3) 3*Y(j,3) 0]);end
out.tetrahedron=struct('family',{T});[~,out.tetrahedron.result] = qstab.zeroexclude(T,'MaxTime',Inf);
[~,out.tetrahedron.witness] = qstab.witnessData(T,out.tetrahedron.result);
VD={[1 -1.2 0],[1 -1.4 0]};out.delay=struct('family',{VD});
[~,out.delay.result] = qstab.zeroexclude(VD,'Variable','zinv','MaxTime',Inf);
[~,out.delay.witness] = qstab.witnessData(VD,out.delay.result);
loss={[1 1],[-1 -1]};out.degreeLoss=struct('family',{loss});
[~,out.degreeLoss.result] = qstab.zeroexclude(loss,'MaxTime',Inf);
[~,out.degreeLoss.witness] = qstab.witnessData(loss,out.degreeLoss.result);
[~,out.lyapunov.result] = qstab.robuststable({[1 1],[1.1 1]},'StoreCertificates',true,'MaxTime',Inf);
[~,out.lyapunov.coverage] = qstab.certificateData(out.lyapunov.result);
out.ok=~out.bernstein.coarse.zeroExcludedNumerically&&out.bernstein.refined.zeroExcludedNumerically&& ...
    strcmp(out.bernstein.result.status,'stable')&&strcmp(out.correlated.result.status,'stable')&& ...
    strcmp(out.interval.result.status,'stable')&&strcmp(out.budget.result.status,'undetermined')&& ...
    out.tetrahedron.witness.isConfirmedCounterexample&&out.delay.witness.isConfirmedCounterexample&& ...
    out.degreeLoss.witness.isDegreeLoss;
if o.Plots
    [~,h] = qstab.plotEnclosure(out.bernstein.coarse,'Samples',out.bernstein.samples,'Visible',o.Visible);figures(end+1)=h.figure;
    [~,h] = qstab.plotEnclosure(out.bernstein.refined,'Samples',out.bernstein.samples,'Visible',o.Visible);figures(end+1)=h.figure;
    [~,h] = qstab.plotEnclosure(out.bernstein.leaf,'Visible',o.Visible);figures(end+1)=h.figure;
    [~,h] = qstab.plotCertificate(out.bernstein.result,'Visible',o.Visible);figures(end+1)=h.figure;
    [~,h] = qstab.plotCertificate(out.budget.result,'Visible',o.Visible);figures(end+1)=h.figure;
    [~,h] = qstab.plotWitness(out.tetrahedron.witness,'Dimension',3,'Visible',o.Visible);figures(end+1)=h.figure;
    [~,h] = qstab.plotWitness(out.delay.witness,'Visible',o.Visible);figures(end+1)=h.figure;
    [~,h] = qstab.plotCertificate(out.lyapunov.result,'Visible',o.Visible);figures(end+1)=h.figure;
end
if o.Verbose
    fprintf('qstab diagnostic examples: ok=%d\n',out.ok);
    fprintf('  Bernstein coarse hull contains zero; refined local hulls exclude zero.\n');
    fprintf('  Stored cells: polynomial %d, interval %d, correlated %d.\n', ...
        out.bernstein.coverage.storedAcceptedCount,out.interval.coverage.storedAcceptedCount,out.correlated.coverage.storedAcceptedCount);
    fprintf('  Partial reduction: %s; unvisited obligations remain visible.\n',out.budget.result.status);
    fprintf('  Tetrahedron witness: %s; delay witness: %s.\n',out.tetrahedron.witness.checkStatus,out.delay.witness.checkStatus);
    fprintf('  Degree loss: %s (not a finite-root instability claim).\n',out.degreeLoss.witness.checkStatus);
end
end
function P=poly(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
