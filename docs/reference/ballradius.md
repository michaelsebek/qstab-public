# qstab.ballradius

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/ballradius.m)

## Help and calling sequences

```text
BALLRADIUS Global scalar weighted-coefficient-ball stability radius.
  [BOUNDS,Info] = qstab.ballradius(P,weights,...)
  BOUNDS = [lower upper] gives the computed global radius bracket;
  it is not the assertion of a single exact radius. Interpret the values
  together with Info.status, witnesses and termination diagnostics.
  Info is the full diagnostic structure returned by earlier releases.
[~,R] = qstab.ballradius(P,GAMMA,'Exponent',ETA,'Variable','s'|'z'|'zinv')
considers full quaternion perturbations delta_k with
      norm((abs(delta_k)/GAMMA(k)),ETA) <= r.
Zero weights fix their coefficients. ETA in [1,Inf] defaults to 2.

The qB02 class-minimum modulus divided by its Holder envelope is minimized
on the stability boundary. CT uses t=w/(1+w) and a homogeneous Bernstein
enclosure on [0,1], INCLUDING the infinite-frequency endpoint. DT uses
the complete upper unit semicircle and a matrix Taylor enclosure. Delay
coefficients and weights are reversed at the DECLARED grade before DT
testing; reconstructed members are returned in the original input variable.

R.lower/R.upper are numerical global bounds even on budget exhaustion.
A grid is not used to assert the lower bound. The upper bound is supported
by a reconstructed finite-boundary member or an explicit degree-loss
member. Converged means gap<=AbsTol+RelTol*upper, NOT validated arithmetic.
CT leading-coefficient loss is not fabricated as a finite-frequency zero.

Options: Exponent=2, Variable='s', Grade=[], AbsTol=1e-5, RelTol=1e-4,
Tolerance=1e-10, MaxCells=8191 (total constructed cells), MaxDepth=40,
MaxTime=30, StoreCells=true. No custom classes or external solver is used.

  Examples (ascending coefficients):
  [bounds,R] = qstab.ballradius([2 1],[1 0],'AbsTol',0.05,'RelTol',0.01);
  disp(bounds);

  Interpretation:
  Weighted full-quaternion coefficient balls, not component boxes. The continuous infinite-frequency endpoint can be degree loss without a finite boundary zero. Inspect lower/upper/limitingEvent and the remaining budget.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `Exponent` | `2` | Block norm exponent eta in [1,Inf]. eta=Inf gives independent full quaternion balls; zero weights fix coefficients. |
| `Variable` | `'s'` | Time convention: s (continuous), z (forward shift) or zinv (delay). Reinterprets the supplied coefficients; it does not transform a system. |
| `Grade` | `[]` | Declared nonnegative integer grade, at least the supplied coefficient layout; [] infers it. Pads exact zeros, never discards coefficients. qB01 Grade applies only to discrete branches. |
| `AbsTol` | `1e-5` | Absolute component of the requested global bracket gap. Converged does not mean an exact radius. |
| `RelTol` | `1e-4` | Relative component of the requested gap, scaled by the finite upper bound. |
| `Tolerance` | `1e-10` | Numerical residual/separation or pairing guard as described for this routine. It is not a physical uncertainty radius and not a directed-rounding proof. |
| `MaxCells` | `8191` | Finite computation budget. Meaning is total constructed/processed covering cells in the selected solver; exhausted cells remain unresolved. |
| `MaxDepth` | `40` | Maximum subdivision depth. It is local to a branch/cell, not a promise of full coverage. |
| `MaxTime` | `30` | Time budget in seconds; Inf disables it when accepted. No incomplete run is relabelled stable. |
| `StoreCells` | `true` | Retain covering cells in diagnostics. Turning storage off does not remove the solver's coverage obligations. |
