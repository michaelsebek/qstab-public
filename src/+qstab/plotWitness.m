function [ax,Info] = plotWitness(W,varargin)
%PLOTWITNESS Draw a previously reconstructed member's algebraic zero classes.
%   [AX,Info] = qstab.plotWitness(W,...)
%   AX is the target MATLAB axes. Info contains all previous plot handles
%   and the plotted data or diagnostics. Request Info for reuse or inspection.
%   Info is the full graphics-handle structure returned by earlier releases.
% [~,H] = qstab.plotWitness(W) uses qstab.witnessData output; it does NOT rerun roots
% or any robust test. Representation='input' (default) or 'forward' explicitly
% selects the delay/forward variable. Other options: Dimension=2|3, Parent,
% Visible. A failed membership check is rejected, never plotted as a witness.
% A degree loss is annotated without manufacturing finite unstable roots.
%
%   Examples (ascending coefficients):
%   F = {[-1 1],[-2 1]}; [~,R] = qstab.robuststable(F);
%   [~,W] = qstab.witnessData(F,R);
%   [ax,H] = qstab.plotWitness(W,'Dimension',2);
%
%   Interpretation:
%   Rendering only: input and forward delay coordinates remain distinct. An unavailable witness is not invented.
%
%   Reference: docs/reference/plotWitness.html.
%   See also qstab.about.
try
Info = qstab.internal.plotWitnessReport(W,varargin{:});
ax = Info.axes;
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
