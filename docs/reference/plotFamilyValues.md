# qstab.plotFamilyValues

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/plotFamilyValues.m)

## Help and calling sequences

```text
PLOTFAMILYVALUES Plot values of sampled family members at one frequency.
  qstab.plotFamilyValues(F,omega)                100 sampled members
  qstab.plotFamilyValues(F,omega,N)              N sampled members
  [ax,D] = qstab.plotFamilyValues(F,omega,N,...) plot and keep the data
  [~,D] = qstab.plotFamilyValues(F,omega,N,'Plot',false)  no figure

  For every sampled P(v)=sum a_k v^k (coefficients on the left), compute
  its RIGHT point evaluation w=sum a_k*xi^k. Choose a pure unit quaternion
  u with 'AxisUnit' (default i). For variable s use xi=u*omega, where
  omega is in radians per unit time. For z use cos(omega)+u*sin(omega),
  for delay use cos(omega)-u*sin(omega); omega is in radians PER SAMPLE
  in both discrete cases. Signed real omega is accepted. To query physical
  frequency Omega in rad/s at sample period h, pass omega=Omega*h.

  'Projection'   'class' (default): (Re w,norm(Im w)). This preserves
                 distance to the origin and maps only zero to the origin;
                 directions in the imaginary 3-space are not preserved.
                 'slice': (Re w,Im(w) dot u), colored by the length of
                 its discarded imaginary component perpendicular to u.
                 A point at the plotted origin may have w ~= 0.
                 'normpoly': planar complex N_P(zeta)=sum_l h_l(zeta)^2
                 for the FOUR real component polynomials h_l, at central
                 complex zeta. This tests algebraic zero classes; it is
                 not the quaternion value w or its squared norm.
  'AxisUnit'     'i','j','k', a real unit 3-vector, or a scalar pure unit
                 MATLAB quaternion or a,b,c,d component struct. Fixes the evaluation direction;
                 the norm-polynomial plane is independent of this choice.
  'AxisUnit','sphere' with 'Projection','class' samples K input directions
                 per polynomial: +i,+j,+k,-i,-j,-k first (as many as fit),
                 then a Fibonacci sphere. 'DirectionCount',K defaults to 32.
                 N counts members; N*K <= 1e6. Values at real arguments
                 repeat across all directions. No sphere search is needed
                 in normpoly mode. A signed slice requires one fixed axis.
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

  'Sampling'     'grid' (default) or 'random'; N counts family members.
  'Seed'         Integer 0..2^32-1 for repeatable random sampling with
                 no change to the caller's RNG. Without Seed, random
                 sampling advances the caller's RNG.
  N is at most 1e6, with at most 1e7 coefficient/parameter samples.
  'ParameterBox' Finite 4-by-2 component limits for a oneparameter model;
                 required for an unrestricted t; a subbox for bounded t.
  'Plot'         true (default) or false; false returns ax=[].
  'Parent'       Existing axes; 'Visible','off' hides a new figure.

  D.valueComponents has N*K rows (K=1 by default), with columns
  [real i j k] of w. D.memberIndex and D.directionIndex identify each row.
  D.valueNorm is the Euclidean norm of each w. D.projected has two
  columns. D.offSlice records the
  perpendicular imaginary length. D.normValues is the complex N-by-1
  array in 'normpoly' mode, empty otherwise. D.complexArgument is zeta,
  D.argumentComponents lists the K query quaternions xi; D.axisUnit lists
  their K pure unit directions as real 3-vectors.
  D.sample gives weights/parameters; D.coefficients(:,k+1,n) contains the
  four real parts of member n at power k. D.isCertificate=false: these
  are samples, not a certified value enclosure or zero-exclusion test.

  Example (ascending coefficients; no add-ons needed):
    minusJ = struct('a',[0 1],'b',[0 0],'c',[-1 0],'d',[0 0]);
    F = qstab.family('polytope',{[0 1],minusJ});
    qstab.plotFamilyValues(F,1,41,'AxisUnit','j');          % actual P(j)
    qstab.plotFamilyValues(F,1,41,'Projection','normpoly'); % central N_P(1i)
  See also qstab.plotFamilyZeros, qstab.valueset, qstab.plotValueSet, qstab.valueenclosure.
```
