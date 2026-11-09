# qstab.plotCertificate

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/plotCertificate.m)

## Help and calling sequences

```text
PLOTCERTIFICATE Plot retained coverage descriptors or their separation gaps.
  [AX,Info] = qstab.plotCertificate(R,...)
  AX is the target MATLAB axes. Info contains all previous plot handles
  and the plotted data or diagnostics. Request Info for reuse or inspection.
  Info is the full graphics-handle structure returned by earlier releases.
[~,H] = qstab.plotCertificate(R,'View','coverage'|'separation'|'auto')
accepts a robust-test report or qstab.certificateData output. No solver
is called. Coverage bars show each ORIGINAL stored record on its own row:
they are projections onto frequency, NOT proof of parameter coverage.
Unvisited subfamily descriptors lacking frequency bounds are explicitly
counted. Missing retained certificates never masquerade as complete data.
Separation values retain their units and are not robustness margins.
Options: Parent, Visible, RowIndices=[], MaxRows=200.

  Examples (ascending coefficients):
  [~,R] = qstab.zeroexclude({[1 1],[2 1]},'StoreCertificates',true);
  [ax,H] = qstab.plotCertificate(R);

  Interpretation:
  Rendering only: unresolved regions and unvisited subfamilies remain visible. Stored evidence, not a fresh optimization, is shown.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `View` | `'auto'` | coverage, separation or auto. Covering bars are frequency projections, not a new proof of parameter coverage. |
| `Parent` | `[]` | Existing axes handle; its current content and hold/NextPlot state are preserved by the public plotting routine. |
| `Visible` | `'on'` | on/off for a newly created figure. Does not alter the data or the validity of a result. |
| `RowIndices` | `[]` | Distinct retained rows to display; empty uses the default subset. Omitted rows remain explicit. |
| `MaxRows` | `200` | Maximum number of displayed certificate rows, not an algorithmic computation budget. |
