# example_qstab_oneparameter_cases

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf)

## Help and calling sequences

```text
QSTAB_ONEPARAMETER_CASES Data for the scalar one-parameter article, revision 7.
  C=qstab_oneparameter_cases returns the example catalogue without running
  a solver. Coefficients in C.data are explicit four-component arrays;
  they are converted to native quaternion inputs by run_qstab_oneparameter.
  C.rays contains the rational directions, unchanged in length or side.
  Recorded analytical targets are references, not newly computed answers.

  E1--E7 are the seven archived frequency datasets. E5V is the additional
  noncommuting tangency, E5C the real cubic multiplicity diagnostic, and
  D11--D13 the three remaining directional diagnostics. A01 checks zero
  perturbations, nonmonic data and nonunit directions.
  See also run_qstab_oneparameter, plot_qstab_oneparameter.
```
