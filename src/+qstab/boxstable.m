function [tf,Info] = boxstable(P,widths,r,varargin)
%BOXSTABLE Strict stability of a componentwise box about a scalar polynomial.
%   [TF,Info] = qstab.boxstable(P,widths,r,...)
%   TF is a double: 1 strictly stable, 0 demonstrably not strictly stable
%   (including a confirmed boundary member), NaN unresolved or prerequisites
%   not met. Inspect Info.status to distinguish these cases; test isnan(TF)
%   before using TF in a branch. Info retains the complete diagnostic report.
%   Info is the full diagnostic structure returned by earlier releases.
% [~,R] = qstab.boxstable(P,W,RADIUS,'Variable','s'|'z'|'zinv') tests ALL real
% coefficient components satisfying abs(delta_knu)<=RADIUS*W_knu. W is a
% real nonnegative 4-by-(grade+1) array (or explicit a,b,c,d component struct).
% Zero components of W are fixed, not independent quaternion balls.
%
% The verified whole-family solver supplies stable decisions. Unstable
% witnesses are reconstructed and spectrally rechecked in the original box.
% Two original-box corners are negative-only candidate checks, never a
% vertex-only certificate. Every moved candidate is spectrally rechecked.
% Padding of floating-point endpoints can only make positive tests harder;
% a witness outside the requested box is NOT accepted as its counterexample.
% Required leading-coefficient loss is assumptionFailed, not a finite root.
%
% Options: Variable='s', Grade=[], Method='zeroexclusion'|'lyapunov',
% MaxCells=2048, MaxDepth=30, MaxTime=20, Tolerance=1e-10,
% StoreCertificates=false. No new class or additional toolbox is required.
%
%   Examples (ascending coefficients):
%   W = zeros(4,2); W(1,1) = 1;
%   [tf,R] = qstab.boxstable([2 1],W,0.5);
%   disp(tf);
%
%   Interpretation:
%   A complete component box is tested. Checked corners can disprove stability but cannot prove a general family stable. Fixed components and original membership must be preserved.
%
%   Reference: docs/reference/boxstable.html.
%   See also qstab.about.
try
Info = qstab.internal.boxstableReport(P,widths,r,varargin{:});
tf = qstab.internal.strictStabilityDecision(Info);
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
