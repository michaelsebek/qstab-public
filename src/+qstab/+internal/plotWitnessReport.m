function H = plotWitnessReport(W,varargin)
%PLOTWITNESSREPORT Internal implementation returning the full report.
ip=inputParser;ip.FunctionName='qstab.plotWitness';
addParameter(ip,'Representation','input',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'Dimension',2,@(x)isnumeric(x)&&isscalar(x)&&ismember(x,[2 3]));
addParameter(ip,'Parent',[]);addParameter(ip,'Visible','on');parse(ip,varargin{:});o=ip.Results;
rep=lower(char(o.Representation));
if ~ismember(rep,{'input','forward'}),error('qstab:plotWitness:Representation','Use input or forward.');end
if ~isstruct(W)||~isscalar(W)||~isfield(W,'schema')||~strcmp(W.schema,'qstab.witness-data.1')|| ...
        ~W.available||~W.isAdmissible
    error('qstab:plotWitness:Witness','Supply an available, numerically admissible witnessData result.');
end
if strcmp(rep,'input'),lambda=W.inputRepresentatives;variable=W.variable;
else,lambda=W.testRepresentatives;variable=W.testVariable;end
if isempty(lambda)
    [fig,ax,restore]=qstab.internal.diagnosticAxes(o.Parent,o.Visible,'qstab witness diagnosis'); %#ok<ASGLU>
    note=text(ax,.5,.5,W.message,'Units','normalized','HorizontalAlignment','center','Interpreter','none');
    H=struct('figure',fig,'axes',ax,'annotation',note,'isCertificate',false);
else
    if ~W.spectrumValid,error('qstab:plotWitness:Spectrum','The stored spectrum is not valid; no class picture was certified.');end
    H=qstab.internal.plotSpectrumReport(lambda,'Variable',variable,'Dimension',o.Dimension, ...
        'Parent',o.Parent,'Visible',o.Visible,'Tolerance',W.tolerance);
end
title(H.axes,{sprintf('Actual family member | %s | variable %s',W.checkStatus,variable), ...
    'Algebraic classes, not an assertion that every class point is an evaluation zero'},'Interpreter','none');
H.witnessData=W;H.representation=rep;H.variable=variable;H.isCertificate=false;
end
