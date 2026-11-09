# run_qstab_papers

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../examples/run_qstab_papers.m)

## Help and calling sequences

```text
RUN_QSTAB_PAPERS Reproduce SELECTED source-linked qB01/qB02/qB03 examples.
E=run_qstab_papers returns ordinary data only: catalogue/reference inputs,
solver outputs, source checks and execution metadata. No figure, directory,
path, RNG, archive or third-party files are modified. It is not a benchmark.

Options: Articles={'qB01','qB02','qB03'}, Cases={} (all selected articles),
UseQuaternion=true, MaxTime=30 PER SOLVER CALL, MaxCells=8192, Verbose=true.
Empty article/case selections, duplicate IDs and unknown names are errors;
selection is validated BEFORE the first solver call. Output order follows
the catalogue, not the selection list. BOUNDED/UNRESOLVED are kept explicit.
E.ok means no reference contradictions/errors, NOT complete reproduction.
E.fullyMatched additionally requires all selected reference checks matched.
To include E in the ONE standard report archive:
  report=test_qstab('PaperExamples',true);
All original numbered figures/examples are NOT reproduced by this suite.
```
