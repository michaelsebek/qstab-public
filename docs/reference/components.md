# qstab.components

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/components.m)

## Help and calling sequences

```text
COMPONENTS Real components of a numeric or native quaternion array.
  [A,B,C,D] = qstab.components(Q) represents Q = A+B*i+C*j+D*k.
  Q may be a MathWorks quaternion array, a real or complex numeric array,
  or an ordinary scalar struct with real numeric fields a,b,c,d of the
  same size. Complex numeric inputs use the 1,i subalgebra. Outputs are
  real double arrays with the same shape as Q (or its component fields).
  No normalization as a rotation is performed. NaN and Inf are rejected.

  Examples (ascending coefficients):
  Q = quaternion(1,2,3,4);
  [a,b,c,d] = qstab.components(Q);
  disp([a b c d]);

  Interpretation:
  The four outputs are real numeric component arrays of exactly the original shape. They are not normalized rotation parameters.

  See also qstab.about.
```
