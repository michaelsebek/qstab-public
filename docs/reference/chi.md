# qstab.chi

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/chi.m)

## Help and calling sequences

```text
CHI Complex adjoint of a quaternion matrix.
  X = qstab.chi(A), for A=Z+W*j with Z=a+i*b and W=c+i*d, returns
               X = [Z,W;-conj(W),conj(Z)].
  A is a numeric/native-quaternion matrix or component struct a,b,c,d.
  This convention preserves products and conjugate transposes. Each
  right quaternion eigenvalue class is represented by a conjugate pair
  in the spectrum of X. No change of the uncertainty set is made.

  Examples (ascending coefficients):
  [X,I] = qstab.chi([-1 0;0 -2]);
  disp(size(X));

  Interpretation:
  The output is an ordinary complex representation, not a quaternionic matrix. Its block convention is [a+1i*b,c+1i*d;-c+1i*d,a-1i*b].

  See also qstab.about.
```
