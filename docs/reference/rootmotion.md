# qstab.rootmotion

Compute the first and second derivatives of a **simple algebraic zero** with respect to the independent real parameters of a scalar polynomial family.

```matlab
[gradient,Hessian,Info] = qstab.rootmotion(F,lambda);
[gradient,Hessian,Info] = qstab.rootmotion(F,lambda,'Parameters',theta);
```

`F` is a `qstab.family` definition of kind `affine`, `interval`, or a bounded `oneparameter` family. It uses ordinary MATLAB data. `lambda` is a finite complex representative of a simple zero of the real norm polynomial at the chosen expansion point. A repeated norm zero is outside this local calculation even if it represents a well-defined quaternionic class.

The default `theta` is the midpoint of the parameter box. Explicit `Parameters` must be a finite real vector inside that box. Its coordinates have exactly the order and units in `qstab.familycomponents(F)`. For an interval, these are actual coefficient components, not normalized positions in their intervals. A polytope must first be written in an explicit affine coordinate system.

The output convention is

\[
\lambda(\theta+h)=\lambda(\theta)+g^{\mathsf T}h+
\frac12 h^{\mathsf T}Hh+O(\lVert h\rVert^3).
\]

`gradient` is a complex column and `Hessian` is a complex symmetric matrix. Transposes in this formula are nonconjugating. At a stationary imaginary-axis contact, the interval article uses

\[
Q=-\frac12\operatorname{Re}H.
\]

The function differentiates the polynomial identity `N(lambda(theta),theta)=0` analytically. It uses neither an eigensolver nor finite differences. For component polynomial values `p_nu`, the norm is the sum of their **complex squares**, which matters at a nonreal argument.

`Info.valid` reports the root residual and simple-root separation checks. If they fail, the primary outputs are `NaN` and `Info.message` explains why. The option `Tolerance` (default `1e-10`, strictly between zero and one) controls these relative guards. Numerical overflow remains unresolved. The result describes local motion; it does not certify the stability of a finite box or the remainder in a Taylor formula.

## A complete native-quaternion example

Consider the affine family

\[
p(s,t)=s+1+t+\mathbf j,\qquad -0.1\le t\le0.1.
\]

Its norm is \(N_p(s,t)=(s+1+t)^2+1\). The upper complex representative is \(\lambda(t)=-1-t+i\), so at the default midpoint the first derivative is -1 and the second derivative is zero.

```matlab
P = quaternion([1 1],[0 0],[1 0],[0 0]);
F = qstab.family('affine',P,{1},[-0.1 0.1]);
[gradient,Hessian,Info] = qstab.rootmotion(F,-1+1i);
gradient                  % -1
Hessian                   % 0
Info.valid                % true
```

The native polynomial coefficient uses Hamilton j; `1i` in the function call is the ordinary complex class representative. Its use does not assert that this representative is a point-evaluation zero of the quaternion polynomial.

Examples I09 and I10 in `run_qstab_intervals` calculate the OOO and EOO Hessians directly from the manuscript coefficients and compare them with the independently stated rational matrices.

See also [family](family.md), [normpoly](normpoly.md), [roots](roots.md), and [the interval study](../interval_article.md).
