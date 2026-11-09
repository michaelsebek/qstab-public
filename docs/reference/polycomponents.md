# qstab.polycomponents

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/polycomponents.m)

## Help and calling sequences

```text
POLYCOMPONENTS Ascending coefficients of a scalar quaternion polynomial.
  C = qstab.polycomponents(P) returns a real 4-by-(grade+1) matrix, with
  columns containing the a,b,c,d parts of coefficients p_0,...,p_grade.
  P is a numeric/native-quaternion coefficient vector, a 1-by-1-by-(g+1)
  coefficient tensor, a vector cell array of scalar coefficients, or a
  component struct a,b,c,d with the same coefficient layout.
  [C,INFO] = ... also records grade, actualDegree and exact zero padding.
  ...,'Grade',G pads to grade G; it never silently removes coefficients.
  All coefficients are ascending. A raw numeric 4-by-N component matrix
  is deliberately not inferred; use the explicit component struct.

  Examples (ascending coefficients):
  [C,I] = qstab.polycomponents([2 3 1]);
  disp(C); disp(I.grade);

  Interpretation:
  The returned 4-by-(g+1) array is an explicitly documented real representation. A raw 4-by-N input is NOT interpreted as quaternion polynomial coefficients.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `Grade` | `[]` | Declared nonnegative integer grade, at least the supplied coefficient layout; [] infers it. Pads exact zeros, never discards coefficients. qB01 Grade applies only to discrete branches. |
