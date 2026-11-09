# qstab.boxradius

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/boxradius.m)

## Help and calling sequences

```text
BOXRADIUS Global brackets for nested scalar coefficient-COMPONENT boxes.
  [BOUNDS,Info] = qstab.boxradius(P,widths,...)
  BOUNDS = [lower upper] gives the computed global radius bracket;
  it is not the assertion of a single exact radius. Interpret the values
  together with Info.status, witnesses and termination diagnostics.
  Info is the full diagnostic structure returned by earlier releases.
[~,R] = qstab.boxradius(P,W,'Variable','s'|'z'|'zinv') considers independent REAL
perturbations abs(delta_knu)<=r*W_knu at the original declared grade.
W is a real nonnegative 4-by-(grade+1) array, not quaternion-ball weights.

LOWER bounds are obtained ONLY from whole-family robust tests, never from
a frequency/direction grid. UPPER bounds have reconstructed finite-boundary
zonotope witnesses, independently checked unstable members, or an explicit
required-leading-coefficient loss. The latter is a regularity event, not
a finite unstable zero. A finite sample cannot prove an infinite margin.

qB02 gives the exact fixed-leading CT zonotope variational formula. This
implementation combines pointwise gauges with existing full-family tests
and monotonic radius bisection. It does NOT claim an exhaustive direct
minimization over quaternion directions. It may stop with a nonzero gap.

Options: Variable='s', Grade=[], Method='zeroexclusion'|'lyapunov',
AbsTol=1e-4, RelTol=1e-3, Tolerance=1e-10, InitialRadius=1,
MaxIterations=32, MaxCells=8192 (GLOBAL across all robust tests),
MaxDepth=30 (per test), MaxTime=30 (global, checked between operations),
MaxRadius=1e6 (exploration cap, NOT a certified upper bound),
BoundaryPoints=[] (additional actual quaternion points ON the boundary),
MaxFacets=4096 per point, StoreCertificates=false.
Default upper seeds: lambda=0 in CT; +/-1 in z and zinv. No normalization
of supplied points/directions occurs. Uncertain trials never shrink bounds.

  Examples (ascending coefficients):
  W = zeros(4,2); W(1,1) = 1;
  [bounds,R] = qstab.boxradius([2 1],W,'AbsTol',0.05,'MaxCells',512);
  disp(bounds);

  Interpretation:
  Global lower bounds come only from complete-family tests of nested boxes. Actual-point zonotope witnesses give upper bounds, not a direction-grid completeness proof. MaxRadius is a search budget, never an upper certificate by itself.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `Variable` | `'s'` | Time convention: s (continuous), z (forward shift) or zinv (delay). Reinterprets the supplied coefficients; it does not transform a system. |
| `Grade` | `[]` | Declared nonnegative integer grade, at least the supplied coefficient layout; [] infers it. Pads exact zeros, never discards coefficients. qB01 Grade applies only to discrete branches. |
| `Method` | `'zeroexclusion'` | Algorithm selection, restricted to the values listed in this function's help. The robuststable default remains Lyapunov subdivision; zeroexclusion is explicit. |
| `AbsTol` | `1e-4` | Absolute component of the requested global bracket gap. Converged does not mean an exact radius. |
| `RelTol` | `1e-3` | Relative component of the requested gap, scaled by the finite upper bound. |
| `Tolerance` | `1e-10` | Numerical residual/separation or pairing guard as described for this routine. It is not a physical uncertainty radius and not a directed-rounding proof. |
| `InitialRadius` | `1` | Initial positive trial radius for the nested-box search; not a certified margin. |
| `MaxIterations` | `32` | Iteration budget for leading-bound search or nested-box radius refinement, as applicable. |
| `MaxCells` | `8192` | Finite computation budget. Meaning is total constructed/processed covering cells in the selected solver; exhausted cells remain unresolved. |
| `MaxDepth` | `30` | Maximum subdivision depth. It is local to a branch/cell, not a promise of full coverage. |
| `MaxTime` | `30` | Time budget in seconds; Inf disables it when accepted. No incomplete run is relabelled stable. |
| `MaxRadius` | `1e6` | Maximum radius to explore, not an upper certificate for the true stability radius. |
| `BoundaryPoints` | `[]` | Optional actual quaternion boundary points used to obtain verified upper candidates. A finite grid never proves a global lower bound. |
| `MaxFacets` | `4096` | Total facet-tuple budget for intrinsic-rank zonotope gauge and reconstruction at each point; incomplete gauge data remain bounded/unresolved. |
| `StoreCertificates` | `false` | Retain accepted certificates for later read-only inspection. It does not strengthen the mathematical assurance level. |
