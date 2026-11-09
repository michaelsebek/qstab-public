function [bounds,Info] = boxradius(P,widths,varargin)
%BOXRADIUS Global brackets for nested scalar coefficient-COMPONENT boxes.
%   [BOUNDS,Info] = qstab.boxradius(P,widths,...)
%   BOUNDS = [lower upper] gives the computed global radius bracket;
%   it is not the assertion of a single exact radius. Interpret the values
%   together with Info.status, witnesses and termination diagnostics.
%   Info is the full diagnostic structure returned by earlier releases.
% [~,R] = qstab.boxradius(P,W,'Variable','s'|'z'|'zinv') considers independent REAL
% perturbations abs(delta_knu)<=r*W_knu at the original declared grade.
% W is a real nonnegative 4-by-(grade+1) array, not quaternion-ball weights.
%
% LOWER bounds are obtained ONLY from whole-family robust tests, never from
% a frequency/direction grid. UPPER bounds have reconstructed finite-boundary
% zonotope witnesses, independently checked unstable members, or an explicit
% required-leading-coefficient loss. The latter is a regularity event, not
% a finite unstable zero. A finite sample cannot prove an infinite margin.
%
% qB02 gives the exact fixed-leading CT zonotope variational formula. This
% implementation combines pointwise gauges with existing full-family tests
% and monotonic radius bisection. It does NOT claim an exhaustive direct
% minimization over quaternion directions. It may stop with a nonzero gap.
%
% Options: Variable='s', Grade=[], Method='zeroexclusion'|'lyapunov',
% AbsTol=1e-4, RelTol=1e-3, Tolerance=1e-10, InitialRadius=1,
% MaxIterations=32, MaxCells=8192 (GLOBAL across all robust tests),
% MaxDepth=30 (per test), MaxTime=30 (global, checked between operations),
% MaxRadius=1e6 (exploration cap, NOT a certified upper bound),
% BoundaryPoints=[] (additional actual quaternion points ON the boundary),
% MaxFacets=4096 per point, StoreCertificates=false.
% Default upper seeds: lambda=0 in CT; +/-1 in z and zinv. No normalization
% of supplied points/directions occurs. Uncertain trials never shrink bounds.
%
%   Examples (ascending coefficients):
%   W = zeros(4,2); W(1,1) = 1;
%   [bounds,R] = qstab.boxradius([2 1],W,'AbsTol',0.05,'MaxCells',512);
%   disp(bounds);
%
%   Interpretation:
%   Global lower bounds come only from complete-family tests of nested boxes. Actual-point zonotope witnesses give upper bounds, not a direction-grid completeness proof. MaxRadius is a search budget, never an upper certificate by itself.
%
%   Reference: docs/reference/boxradius.html.
%   See also qstab.about.
try
Info = qstab.internal.boxradiusReport(P,widths,varargin{:});
bounds = [Info.lower Info.upper];
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
