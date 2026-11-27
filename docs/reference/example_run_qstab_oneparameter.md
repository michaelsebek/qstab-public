# example_run_qstab_oneparameter

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf)

## Help and calling sequences

```text
RUN_QSTAB_ONEPARAMETER Reproduce the scalar one-parameter article examples.
  R=run_qstab_oneparameter runs E1--E7 and the directional diagnostics.
  R=run_qstab_oneparameter('E1') runs the calibration example only.
  R=run_qstab_oneparameter('Cases',{'E1','E5','E5V'}) selects examples.
  C=run_qstab_oneparameter('list') lists all example keys.

  Options (name/value pairs):
    Cases             'all' or a cell of catalogue keys.
    UseQuaternion     true: use MATLAB quaternion; false: components.
    ExactCertificates 'auto' (default), true or false. Auto executes the
                      optional exact certificates when Symbolic Math
                      Toolbox is usable. Otherwise its absence is recorded.
                      No Python process is started by this MATLAB runner.
    Plots             true: save the E1 curves and E3 gain geometry.
    Visible           'off' (default) or 'on'. Existing figures are kept.
    MaxTime           5 seconds per numerical global-radius calculation.
    MaxCells          1024 per global-radius calculation.
    AbsTol, RelTol    1e-4 for numerical global-radius brackets.
    Save              true: MAT/TXT/CSV and an automatic new result ZIP.
    Verbose           true: print example summaries and result locations.
    OutputDir         []: <package>/oneparameter_results/<unique run>.

  The fixed ZIP is <package>/qstab_oneparameter_results.zip; its absolute
  path is R.archiveFile. Errors within individual examples are retained
  in the report and do not prevent the final archive from being written.
  An invalid option or unavailable requested quaternion class is a
  preflight error and creates no report.

  R.ok reports that the executed checks passed. R.articleComplete is a
  stricter field: all examples and all required exact proofs must have run.
  A numerical bound, an unresolved ray contact, or an archived certificate
  is never silently promoted to a newly verified exact result. Timing and
  SVD witness errors are diagnostics, not portable reference constants.
  E7 uses an explicit amplitude/frequency normalization, recorded in R.

  Example:
    R=run_qstab_oneparameter('all','Plots',true);
    R.status
    R.archiveFile
    assert(R.ok,'Inspect the report: %s',R.archiveFile);

  See also qstab_oneparameter_cases, qstab.freqmargin, qstab.raymargin,
  qstab.stabilityradius, qstab.radiuscertificate, qstab.raycertificate.
```
