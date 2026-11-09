function [bounds,Info] = ballradius(P,weights,varargin)
%BALLRADIUS Global scalar weighted-coefficient-ball stability radius.
%   [BOUNDS,Info] = qstab.ballradius(P,weights,...)
%   BOUNDS = [lower upper] gives the computed global radius bracket;
%   it is not the assertion of a single exact radius. Interpret the values
%   together with Info.status, witnesses and termination diagnostics.
%   Info is the full diagnostic structure returned by earlier releases.
% [~,R] = qstab.ballradius(P,GAMMA,'Exponent',ETA,'Variable','s'|'z'|'zinv')
% considers full quaternion perturbations delta_k with
%       norm((abs(delta_k)/GAMMA(k)),ETA) <= r.
% Zero weights fix their coefficients. ETA in [1,Inf] defaults to 2.
%
% The qB02 class-minimum modulus divided by its Holder envelope is minimized
% on the stability boundary. CT uses t=w/(1+w) and a homogeneous Bernstein
% enclosure on [0,1], INCLUDING the infinite-frequency endpoint. DT uses
% the complete upper unit semicircle and a matrix Taylor enclosure. Delay
% coefficients and weights are reversed at the DECLARED grade before DT
% testing; reconstructed members are returned in the original input variable.
%
% R.lower/R.upper are numerical global bounds even on budget exhaustion.
% A grid is not used to assert the lower bound. The upper bound is supported
% by a reconstructed finite-boundary member or an explicit degree-loss
% member. Converged means gap<=AbsTol+RelTol*upper, NOT validated arithmetic.
% CT leading-coefficient loss is not fabricated as a finite-frequency zero.
%
% Options: Exponent=2, Variable='s', Grade=[], AbsTol=1e-5, RelTol=1e-4,
% Tolerance=1e-10, MaxCells=8191 (total constructed cells), MaxDepth=40,
% MaxTime=30, StoreCells=true. No custom classes or external solver is used.
%
%   Examples (ascending coefficients):
%   [bounds,R] = qstab.ballradius([2 1],[1 0],'AbsTol',0.05,'RelTol',0.01);
%   disp(bounds);
%
%   Interpretation:
%   Weighted full-quaternion coefficient balls, not component boxes. The continuous infinite-frequency endpoint can be degree loss without a finite boundary zero. Inspect lower/upper/limitingEvent and the remaining budget.
%
%   Reference: docs/reference/ballradius.html.
%   See also qstab.about.
try
Info = qstab.internal.ballradiusReport(P,weights,varargin{:});
bounds = [Info.lower Info.upper];
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
