function [tf,Info] = ballstable(P,weights,r,varargin)
%BALLSTABLE Strict stability of a scalar weighted coefficient ball.
%   [TF,Info] = qstab.ballstable(P,weights,r,...)
%   TF is a double: 1 strictly stable, 0 demonstrably not strictly stable
%   (including a confirmed boundary member), NaN unresolved or prerequisites
%   not met. Inspect Info.status to distinguish these cases; test isnan(TF)
%   before using TF in a branch. Info retains the complete diagnostic report.
%   Info is the full diagnostic structure returned by earlier releases.
% [~,R] = qstab.ballstable(P,GAMMA,RADIUS,...) uses qstab.ballradius with the same
% options. A radius below the numerical global lower bound is stable; an
% admissible boundary-zero witness strictly inside the requested ball proves
% notStable. A lost required leading coefficient is assumptionFailed, never
% fabricated finite-root instability. A tolerance band/budget is undetermined.
% No general independent-box or correlated family is substituted for a ball.
%
%   Examples (ascending coefficients):
%   [tf,R] = qstab.ballstable([2 1],[1 0],0.5,'AbsTol',0.05);
%   disp(tf);
%
%   Interpretation:
%   A robust decision is about a complete closed coefficient ball. At the numerical boundary the result can remain undetermined; do not coerce an empty isStable flag to a boolean.
%
%   Reference: docs/reference/ballstable.html.
%   See also qstab.about.
try
Info = qstab.internal.ballstableReport(P,weights,r,varargin{:});
tf = qstab.internal.strictStabilityDecision(Info);
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
