# run_qstab_benchmarks

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../examples/run_qstab_benchmarks.m)

## Help and calling sequences

```text
RUN_QSTAB_BENCHMARKS Matched, bounded scalar robustness comparison.
B=run_qstab_benchmarks returns ordinary structs only; it does NOT write
files, rename archives, change paths, seed RNG, or time plot/export work.
To save it in the usual two-entry archive use:
  report=test_qstab('Benchmarks',true);
  % report.benchmarks is inside qstab_public/qstab_test_results.zip.

Methods: lyapunov, direct, four-generator, kharitonov-cube,
kharitonov-faces, quadratic-faces. Unsupported pairs are recorded as
notApplicable, never silently routed to a different method.
Same family and Tolerance/MaxCells/MaxDepth/MaxTime/MaxSubfamilies apply
to all compared calls. Cells from different engines are NOT equal work.
Full-family reductions retain their native GLOBAL subproblem budgets.

Defaults: 3 repetitions, 1 warm-up, MinSampleTime=.1, MaxBatchCalls=32,
MaxCells=512, MaxDepth=32, MaxTime=1, MaxSubfamilies=32, Tolerance=1e-10.
Explicit tic handles are used because the solvers contain their own
tic/toc budget timers. Every timed output is consumed and retained;
witness checking is outside measured time. Method order is rotated.
A short batch is labelled timingAdequate=false, not silently accepted.
No universal faster-method claim or cross-machine comparison is made.

Cases/Methods select cell vectors of exact IDs. List cases using
qstab_benchmark_cases. Warmups and raw timed results are kept in B.rows.
Unknown/budget-limited is valid data; runtime errors or contradictions
set B.ok=false. Completion does NOT imply every family was decided.
```
