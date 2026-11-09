function [ax,Info] = plotCertificate(R,varargin)
%PLOTCERTIFICATE Plot retained coverage descriptors or their separation gaps.
%   [AX,Info] = qstab.plotCertificate(R,...)
%   AX is the target MATLAB axes. Info contains all previous plot handles
%   and the plotted data or diagnostics. Request Info for reuse or inspection.
%   Info is the full graphics-handle structure returned by earlier releases.
% [~,H] = qstab.plotCertificate(R,'View','coverage'|'separation'|'auto')
% accepts a robust-test report or qstab.certificateData output. No solver
% is called. Coverage bars show each ORIGINAL stored record on its own row:
% they are projections onto frequency, NOT proof of parameter coverage.
% Unvisited subfamily descriptors lacking frequency bounds are explicitly
% counted. Missing retained certificates never masquerade as complete data.
% Separation values retain their units and are not robustness margins.
% Options: Parent, Visible, RowIndices=[], MaxRows=200.
%
%   Examples (ascending coefficients):
%   [~,R] = qstab.zeroexclude({[1 1],[2 1]},'StoreCertificates',true);
%   [ax,H] = qstab.plotCertificate(R);
%
%   Interpretation:
%   Rendering only: unresolved regions and unvisited subfamilies remain visible. Stored evidence, not a fresh optimization, is shown.
%
%   Reference: docs/reference/plotCertificate.html.
%   See also qstab.about.
try
Info = qstab.internal.plotCertificateReport(R,varargin{:});
ax = Info.axes;
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
