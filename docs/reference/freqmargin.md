# qstab.freqmargin

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/freqmargin.m)

## Help and calling sequences

```text
FREQMARGIN  Smallest |t| producing a boundary zero at a specified frequency.
  [rho,I]=qstab.freqmargin(F,omega,...) accepts a oneparameter qstab.family.
  These are all-H margins. For finite Bounds use ExtendParameter=true
  only if you intend to extend the law beyond the specified domain.

  rho = qstab.freqmargin(P0,P1,omega)
  [rho,Info] = qstab.freqmargin(P0,P1,omega,'Side','right')
  [rho,Info] = qstab.freqmargin(P0,P1,theta,'Variable','z')

  P0 is a nominal strictly stable scalar polynomial. P1 is its perturbation
  polynomial. Enter numeric/native quaternion coefficients in ascending
  order or use a struct of real component arrays a,b,c,d. In continuous
  time require deg(P0)>=1 and deg(P1)<deg(P0). One constant quaternion t
  multiplies P1 on the chosen side for all coefficients.
  omega is a nonnegative real scalar or vector; each query tests the
  whole imaginary quaternion sphere of radius omega(k). In discrete time
  supply angles theta in [0,pi] radians in the forward z variable.

  rho has the shape of the query. rho(k) is the minimum parameter norm
  for that class, not a global minimum over all frequencies. Inf denotes
  an unattainable contact when the calculation establishes this; check
  Info.valid(k) before interpreting a result. For a finite witness,
  Info.parameterComponents(:,k) and Info.boundaryComponents(:,k) are
  four real components in the order real,i,j,k.

  Options:
    Side         'left' (default): P0+t*P1; 'right': P0+P1*t.
    Variable     's' (default), 'z' (forward), 'zinv' (delay).
    Grade        [] infers the declared grade; discrete polynomials only.
    Tolerance    1e-10, positive and less than one; numerical checks.
    CheckNominal true; false skips the nominal stability check when its
                 result is already known.

  Example:
    [rho,I] = qstab.freqmargin([2 1],1,[0 1 2]);
    disp(rho);  % sqrt(4 + [0 1 2].^2)
    t0 = I.parameterComponents(:,1);

  See also qstab.stabilityradius, qstab.raymargin.
```
