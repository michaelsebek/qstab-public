# qstab.boxfreqmargin

[Function index](../api_index.md) · [MATLAB source](../../src/+qstab/boxfreqmargin.m)

## Purpose

Find the smallest multiplier of an independent coefficient box that produces an algebraic zero on the imaginary quaternion sphere of a specified radius. The uncertainty must affect **one Hamilton component**: the real, `i`, `j` or `k` row, although any number of coefficients in that row may vary independently.

This function minimizes over the entire root class at a fixed frequency. It does not minimize over frequencies. The nominal polynomial need not be stable for this calculation; interpreting the result as a stability reserve requires a stable nominal polynomial and a separate check of degree preservation.

## Calling sequences

```matlab
radius = qstab.boxfreqmargin(P,W,omega)
[radius,Info] = qstab.boxfreqmargin(P,W,omega,'Witness',true)
[radius,Info] = qstab.boxfreqmargin(P,W,omega,'Grade',g)
```

`P` is a scalar polynomial with ascending coefficients. It accepts numeric or native MATLAB `quaternion` coefficients, coefficient cells, and ordinary `a,b,c,d` component structs accepted by `qstab.polycomponents`.

`W` is an explicit nonnegative real array with four rows and `g+1` columns, one width for each Hamilton component of each coefficient. At most one row may contain positive entries. Width zero fixes its coefficient component. All-zero widths are accepted.

The radius multiplies the widths:

\[
|\delta_{k\nu}|\le rW_{\nu,k+1}.
\]

`omega` is a nonempty finite nonnegative real scalar or vector in the polynomial frequency-variable units (rad/s for a physical continuous-time variable). The output retains its row or column shape.

## What is computed

Let \(f_\mu(i\omega)\) denote the four ordinary complex component values of the nominal polynomial, and let \(\nu\) identify the uncertain component. The only possible target values for that component are

\[
\zeta_\pm(\omega)=\pm\sqrt{-\sum_{\mu\ne\nu}f_\mu(i\omega)^2}.
\]

These are complex squares, not absolute squares. The original independent coefficient intervals fill a rectangle at this frequency, with halfwidths

\[
W^e(\omega)=\sum_{k\text{ even}}W_{\nu,k+1}\omega^k,
\qquad
W^o(\omega)=\sum_{k\text{ odd}}W_{\nu,k+1}\omega^k.
\]

For either target, the smallest box multiplier is

\[
r_\pm(\omega)=\max\left\{
\frac{|\operatorname{Re}(\zeta_\pm-f_\nu)|}{W^e(\omega)},
\frac{|\operatorname{Im}(\zeta_\pm-f_\nu)|}{W^o(\omega)}
\right\},\qquad
r(\omega)=\min\{r_+(\omega),r_-(\omega)\}.
\]

When a halfwidth is zero, its displacement must be zero as well. A nonzero frozen displacement cannot be made feasible by increasing the multiplier. Arithmetic ambiguity in this equality returns an unresolved result; the user tolerance is never used to project a nonzero displacement to zero.

The formula is exact in real arithmetic. Its implementation uses ordinary floating point and checks reconstructed original coefficients and a right-evaluation zero. It does not supply outward-rounded enclosures or symbolic algebraic certificates.

## Options

| Option | Default | Meaning |
|---|---|---|
| `Grade` | `[]` | Infer the coefficient grade, or pad exact zero coefficients to the supplied nonnegative integer grade. `W` must have the corresponding number of columns. |
| `Tolerance` | `1e-10` | Positive number below one used for numerical reconstruction checks. It is not a physical radius. |
| `Witness` | `false` | Retain reconstructed original polynomials and actual quaternion right zeros. Reconstruction checks run either way. |
| `Variable` | `'s'` | Continuous time only. Other variables are rejected. |

There is no `Evaluation` option: the class calculation is convention independent, and retained actual-point witnesses use right evaluation. It is also distinct from the multiplication side of a shared one-quaternion parameter.

## Outputs and unresolved cases

- `radius(k)` is the minimum of the two numerically resolved branches. `Inf` with `Info.valid(k)==true` means neither target is attainable by any finite multiplier at that frequency.
- `NaN` with `Info.valid(k)==false` preserves an unresolved or unrepresentable calculation. Underflow is never interpreted as zero width, zero radius or proof of unattainability.
- `Info.branchRadii(:,k)` contains the two target costs. `Info.selectedBranch(k)` identifies the minimizing branch, or an available finite upper candidate when the other branch remains unresolved.
- `Info.lower` and `Info.upper` are guarded floating-point estimates. If one branch is unresolved, the primary radius remains `NaN`; a retained finite witness can still give an upper candidate.
- `Info.frequencyData{k}` retains each branch's status, message, scaled target, even/odd gaps, parity costs, actual boundary point and residual. The nominal component values are divided by the nominal coefficient scale and by the common power scale. Width values are divided by their own coefficient scale and the same power scale. This avoids forming overflowing monomials; the final radius includes the ratio of the coefficient scales.
- With `Witness=true`, `Info.witnesses{k}` contains `polynomial`, `coefficientComponents`, `perturbationComponents`, `boundaryComponents`, `classRepresentative`, `weightedNorm` and `evaluationResidual`. The polynomial is an ordinary component struct accepted by other public functions. The boundary point is an actual right zero; it generally differs from the ordinary complex class representative.

An all-zero nominal polynomial has zero frequency distance. This is not a stability verdict. Uncertain leading coefficients are accepted, but a polynomial may lose its declared degree before any finite-frequency contact; a global reserve must include that separate event.

## Example: one-component slice

The polynomial

\[
p(s)=(s+1)^2+3\mathbf i s+3\mathbf j
\]

is perturbed only in the \(\mathbf i\)-component of its linear coefficient. At \(\omega=1\), the first contact is at multiplier \(3-\sqrt5\).

```matlab
P = struct('a',[1 2 1],'b',[0 3 0], ...
           'c',[3 0 0],'d',[0 0 0]);
W = zeros(4,3);
W(2,2) = 1;

[r,Info] = qstab.boxfreqmargin(P,W,1,'Witness',true);
r
Info.valid
Info.witnesses{1}.boundaryComponents
```

The selected actual root is

\[
-\frac23\mathbf j+\frac{\sqrt5}{3}\mathbf k.
\]

For this particular slice an independent analytic argument excludes all other contact frequencies; the function itself makes no such global claim.

## Example: two independent real feedback gains

Here the constant and linear real components vary independently while the nonreal branches are fixed:

\[
p(s)=s^2+(1+2\mathbf j+2\mathbf k)s+1+2\mathbf i,
\qquad W_{0,0}=W_{0,1}=1/10.
\]

```matlab
P = struct('a',[1 1 1],'b',[2 0 0], ...
           'c',[0 2 0],'d',[0 2 0]);
W = zeros(4,3);
W(1,1:2) = 0.1;

u = roots([1 0 8 -4]);
u = real(u(abs(imag(u)) < 1e-12));
[r,Info] = qstab.boxfreqmargin(P,W,sqrt(u),'Witness',true);
r                              % about 5.14320534829038
Info.witnesses{1}.polynomial
```

The critical real constant and linear coefficients both equal the root \(u\) of \(u^3+8u-4=0\) in \((4/9,1/2)\). The actual right root is

\[
\lambda=-\tfrac12u^2\mathbf i+u\mathbf j-u\mathbf k.
\]

For general boxes involving two or more Hamilton components, use `qstab.boxvalue` at an actual quaternion point or the general box stability routines. Applying this one-component formula to multiple components would change the uncertainty problem and is rejected.
