# qstab.boxvalue

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/boxvalue.m)

## Help and calling sequences

```text
BOXVALUE Component-box distance to a zero at an ACTUAL quaternion point.
[R,I]=qstab.boxvalue(P,W,LAMBDA) minimizes max(abs(delta_knu)/W_knu)
subject to right evaluation of P+delta vanishing at LAMBDA. W is an
explicit NONNEGATIVE REAL 4-by-(grade+1) array, or an a,b,c,d component
struct. Zero entries fix the corresponding real components exactly.
P has the usual ascending scalar coefficient representation.

IMPORTANT: unlike ballvalue, this evaluates an ACTUAL quaternion, NOT a
whole similarity class. Complex LAMBDA lies in the 1,i slice. No arbitrary
direction grid certifies the class minimum. A native quaternion array or
an ordinary a,b,c,d point struct can specify other directions explicitly.

R is a numerical point estimate; I.lower/I.upper retain the bracket.
A bounded pointData status does not assert that the optimum was resolved.
The qB02 zonotope gauge is evaluated in its intrinsic real dimension.
At imaginary points the generators are compressed to at most eight, also
for high-degree polynomials. Primal/dual checks are numerical only.
Options: Evaluation='right'|'left', Grade=[], Tolerance=1e-10,
MaxFacets=4096 (total tuple budget per point, recursive reconstruction
included), Witness=false. Query shapes and coefficient order are retained.

  Examples (ascending coefficients):
  W = zeros(4,2); W(1,1) = 1;
  [r,I] = qstab.boxvalue([2 1],W,0,'Witness',true);
  disp(r); disp(I.valid);

  Interpretation:
  LAMBDA is an ACTUAL quaternion point, not just a class representative. Component boxes are anisotropic. A point witness supplies only an upper candidate for a global radius. Zero width fixes the component exactly.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `Evaluation` | `'right'` | right evaluates sum q_k*lambda^k; left evaluates sum lambda^k*q_k. This is not the side of a shared uncertain parameter. |
| `Grade` | `[]` | Declared nonnegative integer grade, at least the supplied coefficient layout; [] infers it. Pads exact zeros, never discards coefficients. qB01 Grade applies only to discrete branches. |
| `Tolerance` | `1e-10` | Numerical residual/separation or pairing guard as described for this routine. It is not a physical uncertainty radius and not a directed-rounding proof. |
| `MaxFacets` | `4096` | Total facet-tuple budget for intrinsic-rank zonotope gauge and reconstruction at each point; incomplete gauge data remain bounded/unresolved. |
| `Witness` | `false` | Retain reconstructed critical members. Availability and validity flags must be checked before use. |
