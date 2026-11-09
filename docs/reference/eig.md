# qstab.eig

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/eig.m)

## Help and calling sequences

```text
EIG Canonical representatives of right quaternion eigenvalue classes.
  LAMBDA = qstab.eig(A) returns n complex representatives with nonnegative
  imaginary parts for an n-by-n quaternion matrix. Repeated entries are
  retained: the vector has algebraic multiplicity n, not just distinct
  classes. A real eigenvalue of the complex adjoint has doubled
  multiplicity, and its copies are paired explicitly.
  [LAMBDA,INFO] reports the full 2n spectrum, conjugate pairing, numerical
  residuals and diagnostic class groups. Check INFO.isValid before using
  the representatives. Grouping by ClassTolerance is numerical, not an
  exact assertion of multiple eigenvalues.
  Options: 'PairingTolerance' (default 1e-7), 'ClassTolerance' (1e-7).
  'Tolerance',T sets both unless overridden by a subsequent option.

  Examples (ascending coefficients):
  [lambda,I] = qstab.eig(diag([-1+1i,-1+2*1i,-1+3*1i]));
  disp(lambda); disp(I.isValid);

  Interpretation:
  This is the RIGHT spectrum A*v=v*lambda. The second output is diagnostics, NOT eigenvectors. Use eigpairs for V,D. Do not replace the existing qeig function or its two-output convention.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `PairingTolerance` | `1e-7` | Relative tolerance for conjugate pairing of the complex embedded spectrum. Inspect pairing diagnostics and isValid. |
| `ClassTolerance` | `1e-7` | Tolerance for diagnostic grouping of nearby classes; entries and algebraic multiplicities are not deleted. |
| `Tolerance` | `sets both; later named options override` | Numerical residual/separation or pairing guard as described for this routine. It is not a physical uncertainty radius and not a directed-rounding proof. |
