# qstab.plotValueSet

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/plotValueSet.m)

## Help and calling sequences

```text
PLOTVALUESET Plot a scalar norm-polynomial sample cloud and mark zero.
  [AX,Info] = qstab.plotValueSet(S, ...)
  AX is the target MATLAB axes. Info contains all previous plot handles
  and the plotted data or diagnostics. Request Info for reuse or inspection.
  Info is the full graphics-handle structure returned by earlier releases.
  [~,H] = qstab.plotValueSet(S) plots the output of qstab.valueset. Columns
  of S.points correspond to S.lambda; rows are parameter samples.

  Options:
    'Parent'       Existing axes. Their contents and hold state survive.
    'Visible'      'on' (default) or 'off' for a newly created figure.
    'LambdaIndex'  Column index/indices to display; default: all columns.
    'Connect'      false (default). true joins samples in their input
                   order, without interpreting the lines as a boundary.

  No convex hull, enclosure, or robust-stability certificate is inferred
  from the plotted points. In particular, a cloud missing zero does not
  prove that the full value set excludes it. H.isCertificate is false.
  See also qstab.valueset, qstab.plotSpectrum.

  Examples (ascending coefficients):
  [~,S] = qstab.valueset({[1 1],[2 1]},1i,'Weights',[1 0;.5 .5;0 1]);
  [ax,H] = qstab.plotValueSet(S);

  Interpretation:
  Rendering only: discrete sample clouds do not prove zero exclusion and are not silently convexified.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `parent` | `[]` | Existing axes handle; its current content and hold/NextPlot state are preserved by the public plotting routine. |
| `visible` | `'on'` | on/off for a newly created figure. Does not alter the data or the validity of a result. |
| `lambdaindex` | `1:size(S.points,2)` | Indices of the stored lambda columns to display; each selected column is a separate parameter-value set. |
| `connect` | `false` | When true, connect samples in their supplied parameter-row order. Joining samples does not certify interpolation or zero exclusion. |
