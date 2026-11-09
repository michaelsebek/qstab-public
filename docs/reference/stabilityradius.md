# qstab.stabilityradius

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/stabilityradius.m)

## Help and calling sequences

```text
STABILITYRADIUS  Bound the first loss over all directions of one quaternion t.
  [bounds,I]=qstab.stabilityradius(F,...) accepts a oneparameter qstab.family.
  These are all-H margins. For finite Bounds use ExtendParameter=true
  only if you intend to extend the law beyond the specified domain.

  bounds = qstab.stabilityradius(P0,P1)
  [bounds,Info] = qstab.stabilityradius(P0,P1,'Side','right')
  [bounds,Info] = qstab.stabilityradius(P0,P1,'Variable','zinv','Grade',g)

  bounds is the 1-by-2 row [lower upper] for the smallest |t| at which
  strict stability or the required regularity is lost. Inspect Info.status
  before using the endpoints. 'converged' means the numerical gap meets
  the specified tolerance; 'bounded' retains an interval (possibly wide).
  An unresolved calculation does not give an exact margin. An identified
  finite upper witness is described in Info, including its parameter
  components when available; degree loss is distinguished from finite zeros.
  This is a scalar polynomial problem with one shared quaternion parameter.

  P0, P1: numeric/native quaternion coefficients in ascending powers, or
  structs with real arrays a,b,c,d. P0 must be strictly stable. For s,
  deg(P0)>=1 and deg(P1)<deg(P0). A declared discrete grade permits
  varying leading coefficients and is preserved on reversal for zinv.

  Options:
    Side        'left' (default): P0+t*P1; 'right': P0+P1*t.
    Variable    's' (default), 'z' or 'zinv'.
    Grade       [] infer, or discrete declared grade (up to 64).
    AbsTol      1e-6; absolute contribution to the target bracket width.
    RelTol      1e-5; relative contribution (target AbsTol+RelTol*upper).
    Tolerance   1e-10; internal numerical checks.
    MaxCells    4096; search/subdivision budget.
    MaxTime     30 seconds; computational time budget.
    MaxDepth    40; discrete subdivision depth limit.
    StoreCells  true; retain discrete search-cell diagnostics.

  Example: p(s)=s+2+t has its nearest boundary member at t=-2.
    [bounds,I] = qstab.stabilityradius([2 1],1);
    disp(bounds); disp(I.status);  % the bracket contains 2

  These are floating-point bounds, not outward-rounded proof intervals.
  See also qstab.freqmargin, qstab.raymargin, qstab.ballradius.
```
