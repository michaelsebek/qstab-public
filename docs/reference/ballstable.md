# qstab.ballstable

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/ballstable.m)

## Help and calling sequences

```text
BALLSTABLE Strict stability of a scalar weighted coefficient ball.
  [TF,Info] = qstab.ballstable(P,weights,r,...)
  TF is a double: 1 strictly stable, 0 demonstrably not strictly stable
  (including a confirmed boundary member), NaN unresolved or prerequisites
  not met. Inspect Info.status to distinguish these cases; test isnan(TF)
  before using TF in a branch. Info retains the complete diagnostic report.
  Info is the full diagnostic structure returned by earlier releases.
[~,R] = qstab.ballstable(P,GAMMA,RADIUS,...) uses qstab.ballradius with the same
options. A radius below the numerical global lower bound is stable; an
admissible boundary-zero witness strictly inside the requested ball proves
notStable. A lost required leading coefficient is assumptionFailed, never
fabricated finite-root instability. A tolerance band/budget is undetermined.
No general independent-box or correlated family is substituted for a ball.

  Examples (ascending coefficients):
  [tf,R] = qstab.ballstable([2 1],[1 0],0.5,'AbsTol',0.05);
  disp(tf);

  Interpretation:
  A robust decision is about a complete closed coefficient ball. At the numerical boundary the result can remain undetermined; do not coerce an empty isStable flag to a boolean.

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
