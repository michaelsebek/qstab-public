# qstab.certificateEnclosure

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/certificateEnclosure.m)

## Help and calling sequences

```text
CERTIFICATEENCLOSURE Reconstruct ONE retained value-space certificate cell.
  [REGIONS,Info] = qstab.certificateEnclosure(data,R,...)
  REGIONS is the cell array of local outer regions, with polygon and controls.
  Use their union at this argument; Info records the scope and coverage.
  Info is the full diagnostic structure returned by earlier releases.
[~,E] = qstab.certificateEnclosure(F,R,'CellIndex',j) reads R.certificates{j}.
For a polytope, F must be the ORIGINAL family: original vertex weights are
reused, including four-generator reductions. Delay coefficients are reversed
at the recorded grade before evaluating the TEST variable z. This operation
is explicit in E.originalVariable / E.variable.
Affine and even/odd interval leaves already store rectangles; F may then be
[] because their geometry is read from the report. Their source data are not
revalidated. Lyapunov certificates have no norm-value image and are rejected.

Bernstein controls or affine centre rectangles are inflated by their saved
whole-frequency variation radius. A CIRCUMSCRIBED disk polygon is used for
the display, so it is not a falsely shrunken inscribed approximation.
This reads/reconstructs geometry; it never reruns the stability solver and
never elevates a floating-point report to a validated interval proof.

  Examples (ascending coefficients):
  V = {[1 1],[2 1]};
  [~,R] = qstab.zeroexclude(V,'StoreCertificates',true);
  [regions,E] = qstab.certificateEnclosure(V,R,'CellIndex',1);

  Interpretation:
  Reconstructs an outer set for an entire retained parameter/frequency cell. It checks compatibility of the original family. Plotting a cell does not add global coverage.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `CellIndex` | `1` | Select a retained certificate cell. The input family must agree with the saved certificate. |
| `DiskSides` | `64` | Number of sides in the circumscribed rendering polygon for a frequency disk; does not shrink the enclosure. |
