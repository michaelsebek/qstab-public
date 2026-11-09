# run_qstab_paper_extensions

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../examples/run_qstab_paper_extensions.m)

## Help and calling sequences

```text
RUN_QSTAB_PAPER_EXTENSIONS Reproduce SELECTED source-linked qB01/qB02/qB03 examples.
E=run_qstab_paper_extensions returns ordinary data only: catalogue/reference inputs,
solver outputs, source checks and execution metadata. No figure, directory,
path, RNG, archive or third-party files are modified. It is not a benchmark.

Options: Articles={'qB01','qB02','qB03'}, Cases={} (all selected articles),
UseQuaternion=true, Verbose=true. These finite checks have no subdivision budget.
An empty article selection, duplicate IDs and unknown names are errors;
Cases={} selects all cases of the requested articles. Selection is validated
BEFORE the first solver call. Output order follows
the catalogue, not the selection list. BOUNDED/UNRESOLVED are kept explicit.
E.ok means no reference contradictions/errors, NOT complete reproduction.
E.fullyMatched additionally requires all selected reference checks matched.
To include E in the ONE standard report archive:
  report=test_qstab('PaperExtensions',true);
All original numbered figures/examples are NOT reproduced by this suite.
```
