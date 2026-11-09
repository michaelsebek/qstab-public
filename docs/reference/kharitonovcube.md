# qstab.kharitonovcube

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/kharitonovcube.m)

## Help and calling sequences

```text
KHARITONOVCUBE Exact qB02 reduction of a CT quaternion interval polynomial.
  A qstab.family definition supplies its variable, grade and parameter model.
[K,INFO]=qstab.kharitonovcube({L,U}) returns an ordinary affine-family
struct (center,generators,bounds), accepted directly by zeroexclude.
Coefficients are ascending. Parameters 1:4 are even real,i,j,k and 5:8
are odd real,i,j,k, all in [0,1]. Some generators can be identically zero.
K is contained in the original interval family and has EXACTLY its norm
value sets at each imaginary-axis frequency, in exact arithmetic.
Stability of 256 cube VERTICES alone is NOT a robustness test.
No corresponding independent-even/odd claim is made for z or zinv.
Use zeroexclude(...,'Reduction','kharitonov-faces') for pair-free faces.

  Examples (ascending coefficients):
  [K,I] = qstab.kharitonovcube({[1 1],[2 1]});
  [tf,R] = qstab.zeroexclude(K,'Uncertainty','affine');
  disp(R.status);

  Interpretation:
  This is the CT eight-parameter continuous cube, not a list of 256 independently sufficient vertex polynomials. Its independent even/odd construction does not apply to a unit-circle boundary.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `Variable` | `'s'` | Time convention: s (continuous), z (forward shift) or zinv (delay). Reinterprets the supplied coefficients; it does not transform a system. |
