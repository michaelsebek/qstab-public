# run_qstab_public

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../run_qstab_public.m)

## Help and calling sequences

```text
RUN_QSTAB_PUBLIC One standalone acceptance run, without an external qeig.
  REPORT=run_qstab_public temporarily restricts the MATLAB search path to
  this package and existing directories under MATLABROOT, runs installation
  checks, the full regression suite and analytic eigenpair references, then
  restores the caller's path and current folder (also on an exception).

  Run from a freshly restarted PUBLIC MATLAB session. This routine does not
  clear loaded classes/functions, reset persistent/global state or remove
  user data. Path restriction is not a security/process sandbox.

  Options (exact names, logical scalars):
    Examples        true: all existing optional correctness/example groups
                    except the long timing benchmark and figure export
    Figures         false: additionally regenerate/export the eight figures
                    and embed the PNG/PDF bytes in report.publicAcceptance
    DependencyScan  true: MATLAB static dependencies of public src entries
    Verbose         true: print ordinary regression/example progress

  No LegacyFunction/EigenOptions override is accepted. qeig is never requested.
  Timing is false; CheckPrimaryVectors and UseQuaternion are true. Native
  quaternion and graphics regressions are required for acceptance.

  The ONLY result archive remains <package root>/qstab_test_results.zip,
  containing report.mat and report.txt. Installation/session evidence is
  report.publicAcceptance inside MAT. An early preflight/runner failure also
  writes a clearly failed report (zero executed tests is never success).
  run_qstab_public does not edit startup.m, userpath or saved MATLAB path.

  The pre-existing test_qstab remains available for targeted/external work.
  This launcher adds no numerical method and changes no solver tolerance.
```
