# export_qstab_paper_figures

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../examples/export_qstab_paper_figures.m)

## Help and calling sequences

```text
EXPORT_QSTAB_PAPER_FIGURES Export SAVED data without recomputing solvers.
X=export_qstab_paper_figures(F) writes PDF/PNG/FIG/CSV and a text manifest
into <package root>/qstab_paper_figures, independently of the current folder.
Options: OutputDir, Formats={'pdf','png','fig','csv'}, Resolution=180,
Overwrite=true. Relative OutputDir is resolved against this package root.
Only the selected figure files and figures_manifest.txt are replaced.
All figures are generated in a staging folder before committing any target.
No test reports/ZIPs are read, written or renamed. No figure handles are
stored in the returned metadata. CSV z-coordinates are zero for 2-D plots.
```
