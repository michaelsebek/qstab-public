# qstab_guide_example

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../examples/qstab_guide_example.m)

## Help and calling sequences

```text
QSTAB_GUIDE_EXAMPLE Run the small, checked examples printed in the public guide.
  qstab_guide_example() or qstab_guide_example('list') returns the catalogue.
  R=qstab_guide_example('G06') runs one example. R.ok means the example's
  assertions passed; a robust stability verdict, if present, is separate.
  A=qstab_guide_example('all') returns a cell array of completed results or
  raises the first failed assertion. No figure, file, path change, download or
  caller variable is created. A fixed local seed is used in the sampling lab.
  Each case below is the exact source of its corresponding guide listing.
```
