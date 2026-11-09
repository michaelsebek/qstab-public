# qstab.lyapoperator

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/lyapoperator.m)

## Help and calling sequences

```text
LYAPOPERATOR Real matrix of X -> A*X + X*A' on quaternion Hermitian X.
  L = qstab.lyapoperator(A) has size (2*n^2-n)-by-(2*n^2-n).
  A is a square native quaternion/numeric array, or a component struct
  with real fields a,b,c,d. This is NOT a quaternion Kronecker sum.

  Coordinates use a Frobenius-orthonormal quaternion Hermitian basis:
  real diagonal entries first; then sqrt(2)*[a,b,c,d] of each upper
  triangular entry, in row order. INFO documents this coordinate order.
  The construction uses the complex adjoint only for multiplication.

  Examples (ascending coefficients):
  [L,I] = qstab.lyapoperator(-eye(2));
  disp(L);

  Interpretation:
  The output is a REAL operator on quaternion Hermitian matrices, not a quaternionic Kronecker sum. Its eigenvalue interpretation is part of the directional test; do not substitute a naive noncommutative Kronecker identity.

  See also qstab.about.
```
