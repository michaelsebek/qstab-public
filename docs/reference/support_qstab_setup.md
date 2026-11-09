# qstab_setup

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../qstab_setup.m)

## Help and calling sequences

```text
QSTAB_SETUP Connect the public namespace and optional examples for this session.
  INFO=QSTAB_SETUP() adds this package root, src, and examples to path.
  INFO=QSTAB_SETUP('Examples',false) does not add the examples directory.
  QSTAB_SETUP without an output also prints the version and fixed report ZIP.

  Public qstab contains no global q-prefixed wrappers. These belong only to
  the separate private qfiles integration. For compatibility, Wrappers=false
  is accepted; Wrappers=true raises qstab:setup:PrivateWrappers BEFORE any
  path changes. No private toolbox is loaded or initialized.

  Only the parent src of +qstab is added. Never add +qstab directly, and do
  not use genpath on the whole distribution. This function does not modify
  startup.m, savepath, global variables, existing qfiles, or class definitions.
  A stale qwrappers path inside THIS installation is disconnected but its
  files are not deleted; the installation checker reports leftover sources.

  Example:
      % Start in the extracted qstab_public directory.
      info=qstab_setup('Examples',true); % silent setup
      [~,R] = qstab.isstable([2 3 1]);
      disp(R.status);

  Result archive: <this root>/qstab_test_results.zip, independent of pwd.
  See also run_qstab_public, test_qstab, qstab_check_installation, qstab.about.
```
