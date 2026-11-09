# Public function and code key

[Guide (PDF)](qstab_manual.pdf) · [Quick start](../README.md) · [Dependency table](installation.md) · [Companion paper](papers.md)

Use qualified names, e.g. `help qstab.raymargin`. These entries show the source help for this release. The guide keys G01–G11 execute the actual code printed in the guide.

## Input and polynomial representations

| Public call | Purpose |
|---|---|
| [qstab.about](reference/about.md) | Version, scope and dependency information for qstab. |
| [qstab.components](reference/components.md) | Real components of a numeric or native quaternion array. |
| [qstab.polycomponents](reference/polycomponents.md) | Ascending coefficients of a scalar quaternion polynomial. |
| [qstab.chi](reference/chi.md) | Complex adjoint of a quaternion matrix. |
| [qstab.companion](reference/companion.md) | Complex companion matrix or generalized companion pencil. |
| [qstab.normpoly](reference/normpoly.md) | Real scalar norm polynomial, with ascending coefficients. |

## Nominal stability, right spectra and algebraic zero classes

| Public call | Purpose |
|---|---|
| [qstab.isstable](reference/isstable.md) | Decide strict nominal stability of a polynomial or state matrix. |
| [qstab.eig](reference/eig.md) | Canonical representatives of right quaternion eigenvalue classes. |
| [qstab.eigpairs](reference/eigpairs.md) | Right quaternion eigenvectors, diagonal representatives and checks. |
| [qstab.roots](reference/roots.md) | Algebraic zero classes of a quaternion polynomial or matrix polynomial. |

## Families, leading coefficients and robust decisions

| Public call | Purpose |
|---|---|
| [qstab.family](reference/family.md) | Define a scalar quaternion polynomial family using ordinary MATLAB data. |
| [qstab.randfamily](reference/randfamily.md) | Draw a repeatable, well-defined scalar quaternion polynomial family. |
| [qstab.familycomponents](reference/familycomponents.md) | Real component description of a scalar polynomial family. |
| [qstab.leadingbound](reference/leadingbound.md) | Nonvanishing of the FULL leading-coefficient uncertainty set. |
| [qstab.kharitonovcube](reference/kharitonovcube.md) | Exact qB02 reduction of a CT quaternion interval polynomial. |
| [qstab.robuststable](reference/robuststable.md) | Sufficient robust stability test for scalar quaternion families. |
| [qstab.zeroexclude](reference/zeroexclude.md) | Adaptive norm-value zero exclusion for scalar polynomial families. |

## One-quaternion margins and matrix operators

| Public call | Purpose |
|---|---|
| [qstab.raymargin](reference/raymargin.md) | First loss of strict stability along one quaternion direction. |
| [qstab.freqmargin](reference/freqmargin.md) | Smallest /t/ producing a boundary zero at a specified frequency. |
| [qstab.stabilityradius](reference/stabilityradius.md) | Bound the first loss over all directions of one quaternion t. |
| [qstab.lyapoperator](reference/lyapoperator.md) | Real matrix of X -> A*X + X*A' on quaternion Hermitian X. |

## Coefficient balls and component boxes

| Public call | Purpose |
|---|---|
| [qstab.ballradius](reference/ballradius.md) | Global scalar weighted-coefficient-ball stability radius. |
| [qstab.ballset](reference/ballset.md) | Pointwise spectral-set membership for a weighted coefficient ball. |
| [qstab.ballstable](reference/ballstable.md) | Strict stability of a scalar weighted coefficient ball. |
| [qstab.ballvalue](reference/ballvalue.md) | Exact-formula pointwise coefficient-ball distance to a zero class. |
| [qstab.boxradius](reference/boxradius.md) | Global brackets for nested scalar coefficient-COMPONENT boxes. |
| [qstab.boxstable](reference/boxstable.md) | Strict stability of a componentwise box about a scalar polynomial. |
| [qstab.boxvalue](reference/boxvalue.md) | Component-box distance to a zero at an ACTUAL quaternion point. |

## Value sets, evidence and plots

| Public call | Purpose |
|---|---|
| [qstab.valueset](reference/valueset.md) | Evaluate scalar quaternion norm-polynomial value sets. |
| [qstab.valueenclosure](reference/valueenclosure.md) | Norm-value OUTER enclosures at one complex argument. |
| [qstab.certificateData](reference/certificateData.md) | Read retained robust-test cells, without solving again. |
| [qstab.certificateEnclosure](reference/certificateEnclosure.md) | Reconstruct ONE retained value-space certificate cell. |
| [qstab.witnessData](reference/witnessData.md) | Reconstruct and independently check a returned family member. |
| [qstab.plotFamilyZeros](reference/plotFamilyZeros.md) | Plot sampled algebraic zero classes of a quaternion polynomial family. |
| [qstab.plotFamilyValues](reference/plotFamilyValues.md) | Plot values of sampled family members at one frequency. |
| [qstab.plotValueSet](reference/plotValueSet.md) | Plot a scalar norm-polynomial sample cloud and mark zero. |
| [qstab.plotSpectrum](reference/plotSpectrum.md) | Plot right-spectral or algebraic quaternion similarity classes. |
| [qstab.plotCertificate](reference/plotCertificate.md) | Plot retained coverage descriptors or their separation gaps. |
| [qstab.plotEnclosure](reference/plotEnclosure.md) | Draw separate local outer sets, controls and optional samples. |
| [qstab.plotWitness](reference/plotWitness.md) | Draw a previously reconstructed member's algebraic zero classes. |

## Setup, checks and worked examples

- [qstab_check_installation](reference/support_qstab_check_installation.md): Read-only integrity and MATLAB resolution checks.
- [qstab_setup](reference/support_qstab_setup.md): Connect the public namespace and optional examples for this session.
- [run_qstab_public](reference/support_run_qstab_public.md): One standalone acceptance run, without an external qeig.
- [test_qstab](reference/support_test_qstab.md): Run the complete qstab regression and native integration suite.
- [export_qstab_paper_figures](reference/example_export_qstab_paper_figures.md): Export SAVED data without recomputing solvers.
- [plot_qstab_paper_figure](reference/example_plot_qstab_paper_figure.md): Plot one SAVED item; never calls a numerical solver.
- [qstab_benchmark_assess](reference/example_qstab_benchmark_assess.md): Check semantics and reconstruct a returned witness.
- [qstab_benchmark_cases](reference/example_qstab_benchmark_cases.md): Deterministic scalar families, not new stability methods.
- [qstab_benchmark_plan](reference/example_qstab_benchmark_plan.md): Explicit applicable/nonapplicable comparison rows.
- [qstab_benchmark_summarize](reference/example_qstab_benchmark_summarize.md): Preserve unresolved and failed observations.
- [qstab_demo](reference/example_qstab_demo.md): Small examples of quaternion nominal and robust stability.
- [qstab_demo_plots](reference/example_qstab_demo_plots.md): Similarity-class plots and sampled norm-polynomial values.
- [qstab_eigen_assess](reference/example_qstab_eigen_assess.md): Independent reference/multiplicity and eigenpair checks.
- [qstab_eigen_cases](reference/example_qstab_eigen_cases.md): Deterministic matrices with known right-class multisets.
- [qstab_eigen_match](reference/example_qstab_eigen_match.md): Bottleneck bijection of canonical complex multisets.
- [qstab_eigen_timing](reference/example_qstab_eigen_timing.md): Matched ONE-output batches, with retained call results.
- [qstab_guide_example](reference/example_qstab_guide_example.md): Run the small, checked examples printed in the public guide.
- [qstab_paper_cases](reference/example_qstab_paper_cases.md): Selected source-linked examples; no solver calls or files.
- [qstab_paper_check](reference/example_qstab_paper_check.md): Assess a source reference without turning unknown into true.
- [qstab_paper_evaluate](reference/example_qstab_paper_evaluate.md): Execute one catalogue item without writes or plotting.
- [qstab_paper_extension_cases](reference/example_qstab_paper_extension_cases.md): Sixteen additional source-linked finite checks.
- [qstab_paper_extension_evaluate](reference/example_qstab_paper_extension_evaluate.md): Additional finite identities and witnesses.
- [qstab_paper_figures](reference/example_qstab_paper_figures.md): Ordinary data for eight selected article illustrations.
- [qstab_validate_paper_figure](reference/example_qstab_validate_paper_figure.md): Validate ordinary plotting data before side effects.
- [run_qb01_discrete](reference/example_run_qb01_discrete.md): Selected analytic discrete one-quaternion margin examples.
- [run_qb01_examples](reference/example_run_qb01_examples.md): Selected exact examples from the corrected qB01 notes.
- [run_qb02_balls](reference/example_run_qb02_balls.md): Selected weighted-coefficient-ball examples, ordinary data only.
- [run_qb02_boxes](reference/example_run_qb02_boxes.md): Component boxes, not coefficient balls or vertex-only tests.
- [run_qb02_extensions](reference/example_run_qb02_extensions.md): qB02 reductions and correlated CT/DT zero-exclusion.
- [run_qb02_qb03_examples](reference/example_run_qb02_qb03_examples.md): Compact reproducible cases for qB02/qB03 methods.
- [run_qb03_four_generator_examples](reference/example_run_qb03_four_generator_examples.md): Separate article-oriented qB03 demo.
- [run_qstab_benchmarks](reference/example_run_qstab_benchmarks.md): Matched, bounded scalar robustness comparison.
- [run_qstab_diagnostics](reference/example_run_qstab_diagnostics.md): Reproducible enclosure, coverage and witness examples.
- [run_qstab_eigen_comparison](reference/example_run_qstab_eigen_comparison.md): Reference-backed, read-only right-spectrum study.
- [run_qstab_paper_extensions](reference/example_run_qstab_paper_extensions.md): Reproduce SELECTED source-linked qB01/qB02/qB03 examples.
- [run_qstab_papers](reference/example_run_qstab_papers.md): Reproduce SELECTED source-linked qB01/qB02/qB03 examples.
