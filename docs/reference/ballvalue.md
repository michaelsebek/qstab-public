# qstab.ballvalue

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/ballvalue.m)

## Help and calling sequences

```text
BALLVALUE Exact-formula pointwise coefficient-ball distance to a zero class.
[R,INFO]=qstab.ballvalue(P,GAMMA,Z,'Exponent',ETA) evaluates the smallest
|| (|delta_k|/gamma_k) ||_eta that permits an algebraic zero in each class
represented by complex Z. Zero weights fix their coefficients. ETA is in
[1,Inf] (default 2). Coefficients/weights are ascending at the same grade.
Z can be a finite numeric array of any shape; R preserves that shape.

This is the qB02 weighted-coefficient-ball spectral-set formula, evaluated
via the smaller singular value of a 2-by-2 complex polynomial value. It is
NOT abs(p(i*omega)) on one arbitrarily chosen quaternion slice. Its small
singular values use SVD, not a subtraction of two nearly equal squares.
'Witness',true reconstructs an admissible critical polynomial per query.
'Grade',G retains zero padding; weights must include all G+1 coefficients.
Z is always in the INPUT polynomial variable. No delay reversal occurs.
INFO.lower/upper are guarded floating-point pointwise estimates, not
outward-rounded interval proofs. Inspect INFO.valid before using results.

  Examples (ascending coefficients):
  [r,I] = qstab.ballvalue([2 1],[1 0],0,'Witness',true);
  disp(r); disp(I.valid);

  Interpretation:
  Z describes a complete similarity class. Full quaternion coefficient balls are isotropic. Witness polynomials may have nonreal coefficients even for a real nominal polynomial. Pointwise validity is not a global stability certificate.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `Exponent` | `2` | Block norm exponent eta in [1,Inf]. eta=Inf gives independent full quaternion balls; zero weights fix coefficients. |
| `Grade` | `[]` | Declared nonnegative integer grade, at least the supplied coefficient layout; [] infers it. Pads exact zeros, never discards coefficients. qB01 Grade applies only to discrete branches. |
| `Witness` | `false` | Retain reconstructed critical members. Availability and validity flags must be checked before use. |
| `Tolerance` | `1e-10` | Numerical residual/separation or pairing guard as described for this routine. It is not a physical uncertainty radius and not a directed-rounding proof. |
