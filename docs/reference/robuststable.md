# qstab.robuststable

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/robuststable.m)

## Help and calling sequences

```text
ROBUSTSTABLE Sufficient robust stability test for scalar quaternion families.
  A qstab.family definition supplies its variable, grade and parameter model.
  [TF,Info] = qstab.robuststable(data,...)
  TF is a double: 1 strictly stable, 0 demonstrably not strictly stable
  (including a confirmed boundary member), NaN unresolved or prerequisites
  not met. Inspect Info.status to distinguish these cases; test isnan(TF)
  before using TF in a branch. Info retains the complete diagnostic report.
  Info is the full diagnostic structure returned by earlier releases.
  [~,R] = qstab.robuststable(DATA,'Uncertainty',KIND,'Variable',VAR) accepts
  the polytope, interval, affine and oneparameter inputs documented by
  qstab.familycomponents. Coefficients are ASCENDING and may be native
  MATLAB quaternion arrays or native component structs with a,b,c,d fields.

  VAR is 's' (Hurwitz), 'z' (Schur), or 'zinv' (delay z^-1). Delay input
  is reversed at its declared degree and tested in z. 'z^-1' and 'delay'
  are aliases. A nonzero leading coefficient is required throughout the
  family in the variable actually tested; for delay this is the original
  constant coefficient. A witnessed loss of this assumption is reported
  separately and is not silently interpreted as a finite unstable root.

  Method='lyapunov' (default) preserves the previous sufficient
  coefficient-ball Lyapunov test. Method='zeroexclusion' dispatches to
  qstab.zeroexclude for scalar intervals, polytopes, correlated affine
  boxes and one-quaternion families in s, z and zinv. A full-family
  degree certificate is required; separate memberwise normalization is
  not performed by the zero-exclusion route.
  With Method='zeroexclusion', Reduction='four-generator' tests complete
  qB03 generator subfamilies. CT Cartesian intervals additionally support
  kharitonov-cube, kharitonov-faces, and quadratic-faces (monic degree 2).
  MaxCells/MaxTime and MaxSubfamilies are global reduction budgets.
  Reduction='none' (default) preserves the direct subdivision choice.

  Options:
    MaxCells          2048  maximum number of processed subdivision cells
    MaxDepth            30  maximum subdivision depth
    MaxTime             20  seconds (checked between cells; Inf allowed)
    Tolerance        1e-10  numerical safety factor, not a spectral margin
    StoreCertificates false retain each accepted leaf and Lyapunov matrix
    Type       'polynomial' only scalar polynomials in this routine

  R.status is 'stable', 'notStable', 'undetermined', or 'assumptionFailed'.
  R.isStable is respectively true, false, [], or []. An unstable/boundary
  point is a counterexample. Exhaustion of a budget or failure to find a
  common Lyapunov certificate is NOT a counterexample.
  R.unresolvedRegions retains the bounds or original vertex-weight
  matrices of unresolved cells, with a reason and depth, independently of
  StoreCertificates. These descriptors refer to the original input family.

  METHOD. For each parameter cell, enclose every coefficient by a ball
      |a_k-a_k^c| <= r_k.
  If l=|a_n^c|>r_n, normalization b_k=a_n^{-1}a_k satisfies
      |b_k-b_k^c| <= (r_k+r_n*|b_k^c|)/(l-r_n).
  The interleaved complex companion has F=F_c+E*D, where E selects its
  last two rows and ||D||_2 <= delta, the Euclidean norm of these bounds.
  A Hermitian X>0 is computed by a base-MATLAB complex-Schur Lyapunov solve.
  With r the norm of the checked equation residual, sufficient conditions
  for the ENTIRE CELL are
      s: 1-r > 2*||X*E||_2*delta,
      z: 1-r > 2*||F_c'*X*E||_2*delta
                   + ||E'*X*E||_2*delta^2.
  These follow directly by bounding the perturbation in the Lyapunov
  inequality. No test of vertices alone is promoted to robust stability.

  Affine real parameter boxes are bisected. A convex hull is split into
  two hulls by replacing either endpoint of a longest vertex pair with
  its midpoint. Their union covers the original hull, including when the
  listed vertices are affinely dependent. Correlations are retained in
  each cell before taking its sufficient norm enclosure.

  ASSURANCE. The inequalities are sufficient in exact arithmetic. This
  implementation uses double precision with residual/safety checks; it is
  NOT validated interval arithmetic or an exhaustive realization of the
  qBarmish face/zero-exclusion theorems. A stable result requires that all
  covering leaf cells pass the sufficient test. Finite-budget failure is
  reported as undetermined. Nominal qstab.isstable diagnoses candidates
  only after a cell cannot be certified; sampling never proves stability.

  No Control System Toolbox, Optimization Toolbox, custom quaternion class,
  or external quaternion toolbox is required by this implementation.

  Examples (ascending coefficients):
  F = {[1 1],[2 1]};
  [tf,R] = qstab.robuststable(F,'Uncertainty','polytope');
  disp(tf);

  Interpretation:
  This tests scalar polynomial families. The default is Lyapunov subdivision. Stable vertices alone do not certify a general quaternion family. The zeroexclusion method and reductions are explicitly selected.

  See also qstab.about.
```
