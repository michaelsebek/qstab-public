# qstab.intervalfaces

[Function index](../api_index.md) · [MATLAB source](../../src/+qstab/intervalfaces.m) · [Combined stability test](zeroexclude.md)

## Purpose

Construct the complete list of continuous Hurwitz interval test families selected by the degree and parity results of the interval article. This function helps inspect the theorem and reuse its families. It does not test stability.

```matlab
F = qstab.family('interval',L,U);
[families,Info] = qstab.intervalfaces(F);
```

`L` and `U` are ascending scalar polynomial coefficients. Each real, i, j and k component of every coefficient varies independently between its bounds. Coefficients may be native MATLAB `quaternion` arrays, ordinary real/complex arrays, or structs with real fields `a,b,c,d`. As an alternative to the family definition, use:

```matlab
[families,Info] = qstab.intervalfaces({L,U},'Variable','s');
```

## Outputs

`families` is a column cell array. Each `families{k}` is an ordinary affine definition produced in the format of `qstab.family`. It contains `center`, eight `generators`, and an 8-by-2 `bounds` array. All eight cube coordinate names retain their meaning even when a generator is zero. An active free coordinate has bounds `[0 1]`; a fixed coordinate has `[0 0]` or `[1 1]`. Inactive coordinates are fixed at zero to avoid duplicate descriptions.

Each element is directly usable:

```matlab
[tf,Detail] = qstab.zeroexclude(families{1});
```

This call tests the first selected continuous family. It does not decide the other families. To decide the original interval, use the combined call below.

The second output records the construction:

| Field | Meaning |
| --- | --- |
| `totalSubfamilies` | Number of selected continuous families. |
| `parameterDimensionBound` | Largest number of free coordinates among those families. |
| `faceCountsByDimension` | Counts indexed by dimension plus one; entry 1 counts vertices. |
| `descriptors{k}` | Bounds, free and fixed coordinate indices, dimension, parity type and face index for `families{k}`. |
| `selection` | The degree-specific rule that selected the families. |
| `cube` | Kharitonov endpoint pieces, active coordinates and the original leading-coefficient condition. |
| `stabilityTestPerformed` | Always false: enumeration alone makes no stability decision. |

Coordinates 1 through 4 are the even real/i/j/k parts; coordinates 5 through 8 are the corresponding odd parts. A face is **pair-free** when it frees at most one of the even and odd coordinates belonging to the same component. `parityType='EOO'`, for example, means one even and two odd free directions in three distinct components.

## Which families are selected?

| Temporal degree | Rule | Number when all eight cube coordinates are active |
| --- | --- | ---: |
| 1 | Vertices | 256 |
| 2 or 3 | Pair-free two-faces | 1536 |
| 4 | EOO/OOO three-faces and uncovered two-faces | 896 |
| At least 5 | Pair-free three-faces | 1024 |

The full quartic count consists of 512 three-faces and 384 EE two-faces. When some coordinates are inactive, the function keeps any lower-dimensional families that are no longer contained in those three-faces. An EO or OO two-face can then be necessary and is not silently dropped. A constant interval is returned as one complete affine family.

These counts refer to continuous families, not to a finite sample of parameter values. A two-face has two continuously varying parameters. Its four corners alone are generally insufficient.

The reduced coordinates describe correlated endpoint interpolations of several original coefficients of the same parity and component. The cube has the same fixed-frequency norm value sets as the original independent coefficient interval. It is not a replacement of independence by arbitrary coefficient correlations.

## Testing the original interval

```matlab
[tf,Report] = qstab.zeroexclude(F, ...
    'Reduction','interval-faces', ...
    'MaxSubfamilies',1536);
```

or equivalently:

```matlab
[tf,Report] = qstab.robuststable(F, ...
    'Method','zeroexclusion', ...
    'Reduction','interval-faces', ...
    'MaxSubfamilies',1536);
```

The combined test first verifies the degree condition over the original coefficient box. It then tests the selected families with global cell, time and family budgets. `tf` is 1 for a complete numerical sufficient certificate, 0 for a checked admissible counterexample, and NaN for an unresolved result or unmet assumptions. `Report` retains the distinction.

The unchanged default `MaxSubfamilies=1024` is smaller than a full degree-two/three selection. Explicitly raise it to 1536 when complete traversal of that selection is intended. A larger budget permits more work; it does not guarantee that floating-point separation will resolve every family. Unvisited families remain in `Report.unresolvedRegions`.

The legacy reduction names retain their earlier behavior. In particular, `kharitonov-faces` uses the general pair-free three-face selection. The specialized monic Hurwitz `quadratic-faces` construction can be more economical, with at most 192 selected families, and remains available.

## Time convention and assumptions

`Variable` is `'s'` by default. This function applies only to continuous-time Hurwitz intervals. A family declared with `'z'` or `'zinv'` is rejected. The independent even/odd rectangle construction must not be transferred to the unit circle; the original coefficient correlations there are different.

Enumeration is permitted even when `Info.cube.degreeInvariant` is false, so that the proposed faces can be inspected. Such an enumeration is not a certificate. The combined stability test rejects a witnessed loss of the original uniform-degree assumption before reporting any face-based verdict.

The degree/parity selection follows the prepublication interval manuscript's sharp-dimension and quartic-parity results. The numerical tests use the established qstab floating-point enclosures and residual checks. This addition is not an implementation of exact rational quantifier elimination or of the manuscript's exact algebraic decision algorithm.

## A small example

```matlab
F = qstab.family('interval',[1.9 2.9 1],[2.1 3.1 1]);
[families,I] = qstab.intervalfaces(F);
I.totalSubfamilies
I.descriptors{1}

[tf,Report] = qstab.zeroexclude(F,'Reduction','interval-faces');
tf
Report.status
```

Only the real component is uncertain in this example, so the deciding set consists of the four continuous edges of its even/odd rectangle. The polynomial coefficients remain in ascending order.

See also [qstab.kharitonovcube](kharitonovcube.md), [qstab.zeroexclude](zeroexclude.md), [qstab.robuststable](robuststable.md), and [qstab.family](family.md).
