# example_verify_qstab_oneparameter_identities

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf)

## Help and calling sequences

```text
VERIFY_QSTAB_ONEPARAMETER_IDENTITIES Recompute 42 exact article identities.
  R=verify_qstab_oneparameter_identities uses Symbolic Math Toolbox to
  verify the E1 feedback model, calibration margins and the E5 variant.
  R.ok is true only after all 42 checks were freshly executed.
  'Verbose',false suppresses the final count. No Python is launched.
  These identities supplement the generic radiuscertificate/raycertificate
  proofs; they do not replace those global and first-contact certificates.
  See also run_qstab_oneparameter.
```
