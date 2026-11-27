# qstab documentation

[Read the public guide (PDF)](qstab_manual.pdf) · [API reference](api_index.md) · [Installation, dependencies and optional integrations](installation.md) · [Article and exact replay](papers.md) · [Validation](../verification/VALIDATION.md)

**Current release:** 0.20.0. The corrected rc.2 exact branch passed all 32 targeted checks and the complete 13-group/139-comparison article study in native MATLAB, including fresh exact proofs. The final distribution preserves its algorithms and tests. [Scalar one-parameter study: commands, proofs, gain geometry and scaling](oneparameter_article.md).

**Historical 0.19.0 interval article and supplement:** [run all 24 source-linked groups, inspect the comparisons, and replay optional exact certificates](interval_article.md). `run_qstab_intervals('all','Plots',true)` saves numerical results and figures in one current-run archive. The numerical baseline of 0.19.0 was validated as candidate 0.19.0-rc.1 in MATLAB R2024b Update 2 on 20 November 2026: 731 regression checks, all 24 interval groups and 110 comparisons passed. This final packaging is not a fresh native run. The optional MATLAB symbolic replay remains separately unvalidated; see the [validation ledger](../verification/VALIDATION.md).

Run `qstab_guide_example('list')` to see the key and `qstab_guide_example('G06')` for one example. The examples are printed from the same executable MATLAB file. Qualified public entry points are under `qstab.` and return the requested mathematical quantity first, followed by diagnostics where applicable. [Project home](../README.md).
