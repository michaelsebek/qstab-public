function [tf,Info] = isstable(X,varargin)
%ISSTABLE  Decide strict nominal stability of a polynomial or state matrix.
%
%   tf = qstab.isstable(P)
%   [tf,Info] = qstab.isstable(P,'Variable','z')
%   [tf,Info] = qstab.isstable(A,'Type','matrix','Domain','discrete')
%   [tf,Info] = qstab.isstable(Pmat,'Type','polynomialMatrix')
%
%   tf is a scalar double: 1 for strictly stable, 0 for demonstrably not
%   strictly stable (including a confirmed boundary zero), NaN if undecided
%   or required assumptions are not met. Use tf==1, tf==0 and isnan(tf).
%   Info.status distinguishes stable, notStable, undetermined and
%   assumptionFailed; Info.message, Info.method and Info.diagnostics give the
%   reason. To inspect the status, request the SECOND output. Near-boundary
%   floating-point results may remain undecided.
%
%   P: a scalar polynomial in ascending powers with numeric or native
%   quaternion coefficients, or a struct with four real arrays a,b,c,d.
%   A: a square state matrix; specify Type='matrix' explicitly, even for
%   a scalar matrix. Pmat: a square coefficient tensor whose pages hold
%   ascending powers; the nominal polynomial-matrix branch is experimental
%   and requires a nonsingular leading coefficient in the test variable.
%
%   Options:
%     Type       'polynomial' (default), 'matrix', or 'polynomialMatrix'.
%     Variable   's' (default), 'z' (forward shift), 'zinv' (delay).
%     Domain     'continuous' or 'discrete' for state matrices; by default
%                infer from the variable convention (continuous when omitted).
%     Grade      [] infer; a nonnegative integer pads missing high powers.
%                The delay test reverses coefficients at this declared grade.
%     Method     'auto' (default), 'coefficient', 'routh', 'schur',
%                'lyapunov', or 'spectral', as appropriate to the input.
%                Scalar auto uses Routh or Schur on the real norm polynomial;
%                matrix auto uses a Lyapunov/Stein test. 'spectral' explicitly
%                asks for algebraic zeros or right eigenvalues.
%     Tolerance  1e-10; numerical separation and residual tolerance.
%     MaxKroneckerOrder 24; order limit for the direct matrix equation solve.
%
%   Examples:
%     tf = qstab.isstable([2 3 1]);                % 1: (s+1)(s+2)
%     [tf,I] = qstab.isstable([0 1]);              % 0: root at 0
%     disp(I.status)
%     tf = qstab.isstable([0 -1;1 0],'Type','matrix');
%
%   See also qstab.robuststable, qstab.zeroexclude, qstab.roots, qstab.eig.
try
Info = qstab.internal.isstableReport(X,varargin{:});
tf = qstab.internal.strictStabilityDecision(Info);
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
