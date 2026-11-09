# qstab_check_installation

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../qstab_check_installation.m)

## Help and calling sequences

```text
QSTAB_CHECK_INSTALLATION Read-only integrity and MATLAB resolution checks.
  A=qstab_check_installation verifies all installed .m files against the
  distributed MANIFEST.sha256 and checks public/example function paths.
  Missing, modified, unexpected and unhashed code are reported separately.
  It never changes MATLAB path, repairs files, initializes a toolbox or
  calls a solver. In particular it never calls an external qeig.

  Options (exact names):
    PackageRoot      this file's directory by default
    CheckResolution  true; require selected functions in this installation
                     and the native quaternion constructor under matlabroot
    DependencyScan   false; optionally call MATLAB's static dependency tool
                     for src/+qstab public entries (not the external study)

  A.ok is an installation check, not a numerical or publication certificate.
  A.integrity lists source problems; A.resolution lists actual paths.
  Static dependencies cannot prove completeness of dynamic dispatch and
  list only detected installed products. A failed/unavailable requested scan
  is not silently converted to success. SHA-256 is not a signed manifest.
  For temporary isolation and execution use RUN_QSTAB_PUBLIC.

  CheckResolution=false is a file-only inspection of another package root;
  no source file in that other tree is executed. DependencyScan requires
  CheckResolution=true. All outputs are ordinary, serializable data.
```
