function [regions,Info] = valueenclosure(data,lambda,varargin)
%VALUEENCLOSURE Norm-value OUTER enclosures at one complex argument.
%   [REGIONS,Info] = qstab.valueenclosure(data,lambda,...)
%   REGIONS is the cell array of local outer regions, with polygon and controls.
%   Use their union at this argument; Info records the scope and coverage.
%   Info is the full diagnostic structure returned by earlier releases.
% [~,E] = qstab.valueenclosure({P1,...,Pm},lambda) uses ALL pairwise Bernstein
% controls sum_nu p_i,nu(lambda)*p_j,nu(lambda), including diagonal controls.
% [~,E] = qstab.valueenclosure({L,U},lambda,'Uncertainty','interval') retains
% the original real parameters. On the imaginary axis independent intervals
% use the exact even/odd coordinate ranges; elsewhere a joint quadratic
% enclosure is used. No independent recombination of correlated E/O parts.
%
% Optional VertexCells is a cell vector of m-by-k nonnegative matrices:
% each column is a convex combination in the ORIGINAL generator list.
% Optional ParameterBoxes is a cell vector of parameter-by-2 sub-boxes.
% Every local outer set is kept separately. Their union is NOT replaced by
% a single hull. Coverage of USER-SUPPLIED subdivisions is not inferred.
% With no subcells the full parameter domain is enclosed at this argument.
%
% lambda is evaluated in the INPUT polynomial variable, including a delay
% polynomial. No reversal or reciprocal is performed by this function.
% This fixed-argument enclosure is NOT a robust-stability test. It has no
% frequency coverage, stable-reference check, or invariant-degree check.
% No roots, eigenvalues, random samples, or stability solvers are called.
%
% Options: Uncertainty='', VertexCells={}, ParameterBoxes={}, Tolerance=1e-10.
% See also qstab.plotEnclosure, qstab.valueset, qstab.certificateData.
%
%   Examples (ascending coefficients):
%   [regions,E] = qstab.valueenclosure({[1 1],[3 1]},1i);
%   disp(regions);  % an outer set, not an instability verdict
%
%   Interpretation:
%   The outer set includes mixed Bernstein controls or joint correlated component bounds. A zero in the hull may be a false alarm. Multiple local hulls form a union, not one newly convexified family.
%
%   Reference: docs/reference/valueenclosure.html.
%   See also qstab.about.
try
Info = qstab.internal.valueenclosureReport(data,lambda,varargin{:});
regions = Info.regions;
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
