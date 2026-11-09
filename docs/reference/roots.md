# qstab.roots

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/roots.m)

## Help and calling sequences

```text
ROOTS Algebraic zero classes of a quaternion polynomial or matrix polynomial.
  [LAMBDA,INFO] = qstab.roots(P) uses an embedded generalized companion
  pencil. Ascending coefficients are required. Scalar coefficient vectors,
  square coefficient tensors, coefficient cells and a,b,c,d structs are
  accepted. LAMBDA is a complex column with nonnegative imaginary parts;
  multiplicities are retained, including real and infinite classes.
  This computes algebraic classes, not orientation-specific evaluation
  zeros. INFO.fullSpectrum contains the complete embedded spectrum.
  INFO.isValid reports successful conjugate pairing and a numerical
  regularity witness. Singular polynomial matrices are not solved here.
  Options: 'Grade',G; 'Variable','s'|'z'|'zinv' (metadata only);
  'PairingTolerance',T (1e-7); 'ClassTolerance',T (1e-7); 'Tolerance',T
  sets both tolerances. Roots remain in the input variable, also for zinv.

  Examples (ascending coefficients):
  [lambda,I] = qstab.roots([2 3 1],'Variable','s');
  disp(lambda); disp(I.isValid);

  Interpretation:
  These are canonical representatives of algebraic zero classes, with multiplicity. They are not orientation-specific evaluation zeros, not the old qroots contract, and not left eigenvalues. A delay output stays in its input variable.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `Grade` | `[]` | Declared nonnegative integer grade, at least the supplied coefficient layout; [] infers it. Pads exact zeros, never discards coefficients. qB01 Grade applies only to discrete branches. |
| `Variable` | `'s'` | Time convention: s (continuous), z (forward shift) or zinv (delay). Reinterprets the supplied coefficients; it does not transform a system. |
| `PairingTolerance` | `1e-7` | Relative tolerance for conjugate pairing of the complex embedded spectrum. Inspect pairing diagnostics and isValid. |
| `ClassTolerance` | `1e-7` | Tolerance for diagnostic grouping of nearby classes; entries and algebraic multiplicities are not deleted. |
| `Tolerance` | `sets both; later named options override` | Numerical residual/separation or pairing guard as described for this routine. It is not a physical uncertainty radius and not a directed-rounding proof. |
