function out = qstab_demo_plots(varargin)
%QSTAB_DEMO_PLOTS Similarity-class plots and sampled norm-polynomial values.
%   OUT = qstab_demo_plots opens linked-by-data 2-D and 3-D views for the
%   s, z, and z^-1 stability conventions, then plots scalar value samples.
%   OUT = qstab_demo_plots('Visible','off') creates invisible figures.
%   Add the qstab src directory and this examples directory to the path.
%
%   The 3-D views are exact Im_k=0 slices of complete similarity classes.
%   The value-set plot is explicitly a sample plot, not a certificate.

visible = 'on';
if mod(numel(varargin),2) ~= 0
    error('qstab:demoPlots:NameValue','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    if ~strcmpi(varargin{k},'Visible')
        error('qstab:demoPlots:Option','The only option is Visible.');
    end
    visible = varargin{k+1};
end
if isstring(visible) && isscalar(visible), visible = char(visible); end
if ~ischar(visible) || ~any(strcmpi(visible,{'on','off'}))
    error('qstab:demoPlots:Visible','Visible must be ''on'' or ''off''.');
end

q0 = quaternion(0,0,0,0);
q1 = quaternion(1,0,0,0);
qI = quaternion(0,1,0,0);
qJ = quaternion(0,0,1,0);

% Diagonal entries have representatives -0.8+0.45i, 0.25+0.8i,
% and -1.3+0.2i. The j component of the second entry is intentional.
A = quaternion(diag([-.8 .25 -1.3]),diag([.45 0 .2]),diag([0 .8 0]),zeros(3));
lambda = qstab.eig(A);
variables = {'s','z','zinv'};
out = struct();
out.lambda = lambda;
out.spectrum = cell(3,2);
for k = 1:numel(variables)
    f = figure('Color','w','Visible',visible,'Name',['qstab: ' variables{k}], ...
        'Position',[100 100 1100 440]);
    a2 = subplot(1,2,1,'Parent',f);
    a3 = subplot(1,2,2,'Parent',f);
    [~,out.spectrum{k,1}] = qstab.plotSpectrum(lambda,'Variable',variables{k}, ...
        'Dimension',2,'Parent',a2);
    [~,out.spectrum{k,2}] = qstab.plotSpectrum(lambda,'Variable',variables{k}, ...
        'Dimension',3,'Parent',a3);
end

P1 = repmat(q0,1,1,3);
P1(:,:,1) = q1+qI;
P1(:,:,2) = 2*q1;
P1(:,:,3) = q1;
P2 = repmat(q0,1,1,3);
P2(:,:,1) = 2*q1+qJ;
P2(:,:,2) = 3*q1;
P2(:,:,3) = q1;
t = linspace(0,1,61).';
[~,out.valueSamples] = qstab.valueset({P1,P2},1i*[0 .5 1 1.5 2], ...
    'Weights',[1-t t]);
[~,out.valuePlot] = qstab.plotValueSet(out.valueSamples,'Visible',visible);
end
