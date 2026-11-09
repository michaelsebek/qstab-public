# qstab.familycomponents

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/familycomponents.m)

## Help and calling sequences

```text
FAMILYCOMPONENTS Real component description of a scalar polynomial family.
  F = qstab.familycomponents(definition) accepts a finite qstab.family struct.
  F = qstab.familycomponents(DATA,'Uncertainty',KIND) accepts:
    'polytope'    DATA = {P1,...,Pm}; the real convex hull of the vertices.
    'interval'    DATA = {L,U}; componentwise bounds, including all four
                  real components of every quaternion coefficient.
    'affine'      DATA.center, DATA.generators (a cell array), and
                  DATA.bounds (number-of-generators by 2, real).
    'oneparameter' DATA.p0, DATA.p1, DATA.side ('left' or 'right'), and
                  DATA.bounds (4 by 2, for real,i,j,k components of t).
  An omitted KIND is inferred from a struct, or is 'polytope' for a cell.

  All polynomials use ASCENDING coefficients. Each polynomial is accepted
  by qstab.polycomponents; the preferred public type is MATLAB quaternion.
  Native structs with real fields a,b,c,d are also accepted. No custom
  class, quaternion constructor, or external quaternion toolbox is used.

  F is an ordinary struct containing only real double arrays and metadata.
  F.vertices has size 4-by-(degree+1)-by-m for a polytope. Other kinds use
      P(theta) = F.base + sum_j theta(j)*F.generators(:,:,j).
  F.bounds contains the actual real parameters used in this expression.
  For intervals, these parameters are the varying coefficient components
  themselves, identified by F.parameterIndex = [component,coefficient].
  The coefficient index in parameterIndex starts at zero.

  Different vertex/generator lengths are padded with zero HIGH-degree
  coefficients. The resulting declared degree is retained. This function
  does not certify a nonzero leading coefficient over the family.

  The oneparameter form is a general affine family, also when p1 changes
  the leading coefficient. The additional deg(p1)<deg(p0) assumption of a
  particular stability-radius theorem is not imposed by this converter.

  Examples (ascending coefficients):
  F = qstab.familycomponents({[1 1],[2 1]},'Uncertainty','polytope');
  disp(F.kind); disp(F.degree);

  Interpretation:
  This is an inspection/conversion representation with real arrays, not a new polynomial object. Preserve common weights and parameters. The function does not itself prove robust stability.

  See also qstab.about.
```
