# qstab.plotFamilyZeros

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/plotFamilyZeros.m)

## Help and calling sequences

```text
PLOTFAMILYZEROS Plot sampled algebraic zero classes of a quaternion polynomial family.
  qstab.plotFamilyZeros(F)                   sample 100 members of F
  qstab.plotFamilyZeros(F,N)                 sample N members (not N zeros)
  [ax,D] = qstab.plotFamilyZeros(F,N,...)    plot and keep the underlying data
  [~,D] = qstab.plotFamilyZeros(F,N,'Plot',false)  compute without a figure

  For each admissible parameter vector, the function forms one polynomial,
  calls qstab.roots at the declared grade, and plots the numerically
  checked finite representatives Re(q) + 1i*norm(Im(q)).
  They represent similarity classes and retain algebraic multiplicity.
  The coordinates stay in the input variable. Strict stability lies to the left
  of the imaginary axis in s, inside the unit circle in z, and outside
  it in the delay variable z^-1. Infinite representatives are counted
  separately, never drawn at an artificial finite coordinate. Green
  circles: stable side; red crosses: unstable side; amber diamonds: on
  or near the boundary (relative display band 100*eps).

  F may be a struct from qstab.family or qstab.randfamily, or any raw
  input accepted by qstab.familycomponents.
  By default, a cell {P1,...,Pm} denotes a polytope. To sample independent
  component bounds {L,U}, pass 'Uncertainty','interval'. Affine and
  oneparameter definitions are structs (see qstab.familycomponents).
  A canonical qstab.familycomponents result is accepted as well.
  'Variable'      's' (default), 'z' or 'zinv'; matching family metadata
                  are inferred for tagged family structs.
  'Grade'         For raw inputs, pad to this grade (never discard coefficients).
  'Uncertainty'   'polytope', 'interval', 'affine' or 'oneparameter'.

  'Sampling'     'grid' (default) or 'random'. N is the total case budget;
                 grid mode uses a deterministic selection when needed.
  'Seed'         Integer 0..2^32-1 for a repeatable random draw with
                 no change to the caller's RNG. Without Seed, random
                 sampling advances the caller's RNG.
  N is at most 1e6, with at most 1e7 coefficient/parameter samples.
  'ParameterBox' Finite 4-by-2 [lower upper] table of quaternion parameter
                 components, required for an unrestricted oneparameter F.
                 For a bounded oneparameter F this can restrict the domain.
  'Plot'         true (default) or false. If false, ax=[].
  'Parent'       Existing axes. Previous graphics and the hold state remain.
  'Visible'      'on' or 'off' for a newly created figure.

  D.lambda is grade-by-N (one column per member). Inf marks an infinite
  class. D.valid, D.status and D.message flag invalid or unresolved columns;
  never use their candidates as evidence of a zero. D.zeroPolynomial marks
  a member with no defined zero list. D.infiniteCount counts validated
  infinite representatives per member (NaN when invalid). D.sample contains
  the actual parameters/weights; D.coefficients(:,k+1,n) holds real,i,j,k
  parts of the coefficient at power k in member n. D.isCertificate=false.
  A display with no unstable points does not prove robust stability.

  Example (ascending coefficients; no add-ons needed):
    F = qstab.family('interval',[1 1],[3 1]);
    [ax,D] = qstab.plotFamilyZeros(F,21);
    D.lambda                         % points on the real interval [-3,-1]
  See also qstab.plotFamilyValues, qstab.roots, qstab.robuststable, qstab.plotSpectrum.
```
