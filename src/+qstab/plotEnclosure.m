function [ax,Info] = plotEnclosure(E,varargin)
%PLOTENCLOSURE Draw separate local outer sets, controls and optional samples.
%   [AX,Info] = qstab.plotEnclosure(E,...)
%   AX is the target MATLAB axes. Info contains all previous plot handles
%   and the plotted data or diagnostics. Request Info for reuse or inspection.
%   Info is the full graphics-handle structure returned by earlier releases.
% [~,H] = qstab.plotEnclosure(E,'Samples',S), E from valueenclosure or
% certificateEnclosure, S from valueset AT THE SAME complex argument.
% Samples are expected to describe the same original family; their argument
% is checked, but their membership is not re-proved by this plot.
% No new solver calls. Controls are not represented as admissible members;
% local polygons are never merged into a global hull. A zero-containing
% OUTER polygon is not evidence of an actual zero or instability.
% Options: Parent=[], Visible='on', RegionIndices=[], MaxRegions=128,
% ShowControls=true, Samples=[]. Omissions are explicit in H and the title.
%
%   Examples (ascending coefficients):
%   [~,E] = qstab.valueenclosure({[1 1],[3 1]},1i);
%   [ax,H] = qstab.plotEnclosure(E);
%
%   Interpretation:
%   Rendering only: attainable samples, controls and outer regions are distinct objects. Empty groups are legal. A drawing is not a new certificate.
%
%   Reference: docs/reference/plotEnclosure.html.
%   See also qstab.about.
try
Info = qstab.internal.plotEnclosureReport(E,varargin{:});
ax = Info.axes;
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
