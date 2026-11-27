# support_test_qstab_oneparameter

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf)

## Help and calling sequences

```text
TEST_QSTAB_ONEPARAMETER Check the new scalar one-parameter functionality.
  R=test_qstab_oneparameter runs the numerical, article-data and optional
  exact-rational regression groups. It needs no private toolbox classes.
  R=test_qstab_oneparameter('ExactCertificates',true) requires usable
  Symbolic Math Toolbox; false explicitly omits that part; 'auto' (default)
  records it as pending when Symbolic Math Toolbox is unavailable.

  Options: ExactCertificates='auto', Verbose=true.
  R.ok means all executed checks passed and explicit requirements were met.
  R.complete additionally requires all exact checks. Neither field claims
  native quaternion acceptance; run RUN_QSTAB_PUBLIC for that separate gate.
  Source examples with native quaternion are run by RUN_QSTAB_ONEPARAMETER.
  An auxiliary Octave execution is always labelled native-MATLAB pending.

  The current report is saved under <package>/oneparameter_test_results;
  a new verified <package>/qstab_oneparameter_test_results.zip is made even
  if a test fails. Its full path is R.archiveFile.
```
