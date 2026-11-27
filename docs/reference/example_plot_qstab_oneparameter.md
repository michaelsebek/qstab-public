# example_plot_qstab_oneparameter

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf)

## Help and calling sequences

```text
PLOT_QSTAB_ONEPARAMETER Plot computed scalar one-parameter article results.
  FILES=plot_qstab_oneparameter(R,'Directory',folder) saves the E1 frequency
  comparison and E3 nearest-gain geometry when those cases are present in R.
  PNG, PDF, MATLAB FIG and accompanying CSV data are written. FILES gives
  their absolute paths. 'Visible','on' also keeps the new figures open.
  Existing figures and axes are not modified.

  E1 uses the computed curves in R, including the restricted-slice quotient.
  Black dots evaluate the reported reference frequencies; their placement
  is not by itself a new certificate of global minimization.
  E3 distinguishes the analytical sphere from the seven computed witness
  checks. Samples alone do not prove the complete sphere characterization.
  See also run_qstab_oneparameter.
```
