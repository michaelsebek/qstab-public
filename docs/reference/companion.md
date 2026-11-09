# qstab.companion

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/companion.m)

## Help and calling sequences

```text
COMPANION Complex companion matrix or generalized companion pencil.
  [A,B,INFO] = qstab.companion(P) constructs a pencil A-lambda*B for
  P(v)=P0+P1*v+...+Pg*v^g, using the complex adjoint of each coefficient.
  P may be a scalar coefficient vector/tensor, a square matrix coefficient
  tensor, a vector cell array {P0,...,Pg}, or a component struct a,b,c,d.
  Coefficients are ascending and act from the left. The pencil preserves
  infinite eigenvalues when the declared leading coefficient is singular.
  A = qstab.companion(P) instead returns the normalized companion matrix;
  this one-output form requires a numerically nonsingular leading block.
  ...,'Grade',G pads with zero high coefficients without silently trimming.
  ...,'Variable','s'|'z'|'zinv' sets metadata only: no reciprocal transform
  is performed here. All spectral points are in the actual input variable.

  Examples (ascending coefficients):
  [A,B,I] = qstab.companion([2 3 1]);
  e = eig(A,B);  % full complex embedded spectrum, not n quaternion classes
  disp(e);

  Interpretation:
  A and B are COMPLEX embedded matrices of order twice the quaternion companion order. They are not quaternion companion matrices and must not be wrapped as such. One output normalizes the leading block; two outputs preserve the generalized pencil.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `Grade` | `[]` | Declared nonnegative integer grade, at least the supplied coefficient layout; [] infers it. Pads exact zeros, never discards coefficients. qB01 Grade applies only to discrete branches. |
| `Variable` | `'s'` | Time convention: s (continuous), z (forward shift) or zinv (delay). Reinterprets the supplied coefficients; it does not transform a system. |
