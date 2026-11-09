# qstab.ballset

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/ballset.m)

## Help and calling sequences

```text
BALLSET Pointwise spectral-set membership for a weighted coefficient ball.
  [MEMBERSHIP,Info] = qstab.ballset(P,weights,z,r,...)
  MEMBERSHIP is pointwise: 1 inside, 0 outside, NaN unresolved.
  It is not a whole-family stability decision.
  Info is the full diagnostic structure returned by earlier releases.
[~,S] = qstab.ballset(P,GAMMA,Z,R,'Exponent',ETA) classifies the requested complex
points only. Inside means some admissible member has this algebraic class;
outside means the pointwise distance exceeds R with a numerical guard.
Values in the tolerance band remain unresolved. No grid certifies the
continuum between its points or robust stability of the whole family.
Z refers to the INPUT polynomial variable, including a delay variable.
Options are passed to qstab.ballvalue. Witness generation is optional.

  Examples (ascending coefficients):
  [membership,S] = qstab.ballset([2 1],[1 0],[-2 -1 0],0.5);
  disp(membership);

  Interpretation:
  Only the supplied query points are classified. A sampled contour or grid never certifies the whole stability region.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `Exponent` | `2` | Block norm exponent eta in [1,Inf]. eta=Inf gives independent full quaternion balls; zero weights fix coefficients. |
| `Grade` | `[]` | Declared nonnegative integer grade, at least the supplied coefficient layout; [] infers it. Pads exact zeros, never discards coefficients. qB01 Grade applies only to discrete branches. |
| `Witness` | `false` | Retain reconstructed critical members. Availability and validity flags must be checked before use. |
| `Tolerance` | `1e-10` | Numerical residual/separation or pairing guard as described for this routine. It is not a physical uncertainty radius and not a directed-rounding proof. |
