# run_qstab_eigen_comparison

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../examples/run_qstab_eigen_comparison.m)

## Help and calling sequences

```text
RUN_QSTAB_EIGEN_COMPARISON Reference-backed, read-only right-spectrum study.
E=run_qstab_eigen_comparison checks qstab.eig on fourteen deterministic
matrices with analytic reference multisets, preserving algebraic counts.
E=run_qstab_eigen_comparison('LegacyFunction',@qeig,'Timing',true)
additionally invokes the EXISTING function explicitly supplied by the user.
Initialize that function's own environment before calling this runner.
No qeig alias, path addition, global initialization or source edit occurs.

Options: Cases (IDs or 'all'), UseQuaternion (true), Tolerance (1e-7),
LegacyFunction ([]), CheckLegacyVectors (true), Timing (false),
CheckPrimaryVectors (false) explicitly checks qstab.eigpairs as well,
outside one-output timing. Its assessment is retained in primary.eigenpairs.
Repetitions (3), WarmupRuns (1), MinSampleTime (.05), MaxBatchCalls (256),
Verbose (true). With no legacy function, timing comparison is not applicable.
Existing [V,D]=qeig(A) is checked separately from qstab's [lambda,info].
A defective fixture does not require a full eigenvector basis.

E.ok requires the selected checks, including requested legacy checks.
COMPLETED_WITH_DIFFERENCES is a completed experiment, NOT an endorsement of
a differing/failed legacy solver. All source data and failures are retained.
Output is ordinary serializable data. This function itself writes no files.
```
