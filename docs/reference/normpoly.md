# qstab.normpoly

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/normpoly.m)

## Help and calling sequences

```text
NORMPOLY Real scalar norm polynomial, with ascending coefficients.
  N = qstab.normpoly(P) forms N(v)=p(v)*conj(p)(v), where conjugation
  acts on coefficients and v is central. If p=a+b*i+c*j+d*k, then
  N=a^2+b^2+c^2+d^2, using polynomial multiplication, not pointwise
  squaring. N has declared grade twice that of P; no zeros are trimmed.
  [N,INFO] = ... returns degree/grade and input metadata.
  ...,'Grade',G has the same meaning as in qstab.polycomponents.

  Examples (ascending coefficients):
  [n,I] = qstab.normpoly([2 1]);
  disp(n);  % [4 4 1], ascending, not MATLAB roots order

  Interpretation:
  Returns a real ascending coefficient vector of declared grade 2*g. Multiplication takes place in the central polynomial ring. N(z) uses complex squares, not absolute squares; N(lambda) is not generally |p(lambda)|^2.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `Grade` | `[]` | Declared nonnegative integer grade, at least the supplied coefficient layout; [] infers it. Pads exact zeros, never discards coefficients. qB01 Grade applies only to discrete branches. |
