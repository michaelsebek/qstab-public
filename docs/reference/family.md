# qstab.family

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/family.m)

## Help and calling sequences

```text
FAMILY Define a scalar quaternion polynomial family using ordinary MATLAB data.
  F = qstab.family('interval',L,U) uses independent bounds on each real,
  i, j and k part of each coefficient. L must be <= U componentwise.
  F = qstab.family('polytope',{P1,...,Pm}) uses one real nonnegative weight
  per generator, summing to one. The weights are common to all powers.
  F = qstab.family('affine',P0,{G1,...,Gm},B) denotes P0+sum(t_j*Gj);
  each real t_j belongs to [B(j,1),B(j,2)], shared across coefficients.
  F = qstab.family('oneparameter',P0,P1,'Side','left','Bounds',B)
  denotes P0+t*P1, with ONE quaternion t. 'right' denotes P0+P1*t.
  B is a real ordered 4-by-2 box on real(t),i,j,k coordinates. Omit B
  or use [] for all of H: valid for a ray margin; a bounded-domain test
  or plot requires an explicit finite domain.

  All input polynomials have ASCENDING coefficients; use MATLAB quaternion,
  numeric real/complex vectors or component structs a,b,c,d. 'Variable'
  is 's' (default), 'z', or 'zinv'. 'Grade' pads with high-order zero
  coefficients; default = maximum supplied grade. 'SamplingPeriod' is 0
  for 's' and positive for discrete variables (default 1). It is metadata;
  no transformation or frequency scaling is applied.
  F is a scalar STRUCT, not a toolbox class. Fields include kind, variable,
  grade, samplingPeriod and lower/upper, vertices, center/generators/bounds,
  or p0/p1/bounds/side according to kind. No stability test is performed.
  F = qstab.family('polytope',{[1 1],[2 1]}); tf=qstab.robuststable(F);
  F = qstab.family('oneparameter',[2 1],1,'Bounds',[-1 .5;zeros(3,2)]);
  [ab,I]=qstab.raymargin(F,-2); I.parameterDomain.admissibleAlpha
  See also qstab.randfamily, qstab.familycomponents, qstab.robuststable.
```
