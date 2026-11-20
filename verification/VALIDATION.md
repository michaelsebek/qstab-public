# Validation of qstab 0.19.0

## Native MATLAB validation of the 0.19.0 sources

**Public acceptance: PASSED.** The release candidate `0.19.0-rc.1` was run on 20 November 2026 with MATLAB **24.2.0.2773142 (R2024b Update 2)** on Windows **PCWIN64**. All **731 regression checks passed**, with zero failures and zero skips. The MAT report, text report and Command Window transcript contain the same cases in the same order. A separate native interval-study run matched **all 24 example groups and 110 comparisons**.

| Scope | Actual evidence |
|---|---|
| Numerical core | 646 regression checks passed, including the new face reductions, frequency margins, analytic root motion and interval-study fixtures. |
| Native types and graphics | All 62 native quaternion/MATLAB checks and all 23 graphics checks passed. `Figures=false` disables an additional publication-figure review; the ordinary graphics regressions still ran. |
| Source identity | SHA-256 and byte sizes for all 208 installed MATLAB files match the delivered `0.19.0-rc.1` archive. The initial inventory, final inventory and report inventory agree. |
| Installation and dependencies | All 84 checked entry points resolved correctly before and after the run. The initial static dependency scan checked 43 public entries and detected 103 files, with none outside qstab and the MathWorks installation. No external eigensolver was requested. Static analysis does not prove coverage of every execution branch. |
| Caller state | The acceptance runner recorded restoration of the caller's MATLAB path and working folder. |
| Interval article and supplement | A separate native quaternion run returned `NUMERICAL_EXAMPLES_MATCHED`: 24/24 groups and 110/110 comparisons. The saved MAT data, text report and comparison CSV agree. Eight figures were exported as PNG, PDF, FIG and CSV; their saved numerical coordinates were checked. |
| Existing selected paper cases | The older catalogues matched 20 selected cases and 16 extension cases. All 14 analytic right-spectrum/eigenpair fixtures matched. These finite catalogues do not claim full reproduction of every earlier paper. |

The regression phase took 80.850 s, the complete acceptance runner 105.929 s and the separate interval study 23.163 s. These are observations of those runs, not performance benchmarks. Successful execution is not a stability conclusion for every example: intentional counterexamples, degree loss, unresolved near-boundary decisions and budget-limited brackets retain their documented meaning.

**From the tested candidate to 0.19.0.** The final release changes version/descriptive text in `qstab_setup.m`, `run_qstab_public.m`, `src/+qstab/about.m` and `src/+qstab/Contents.m`, and two display labels in `examples/plot_qstab_intervals.m`. Numerical algorithms, tolerances, test cases and the plotted values are unchanged. Documentation and citation metadata are updated. The evidence is attributed to the actual tested candidate; no second native execution after the final text edits is claimed.

[Machine-readable native evidence](native-matlab-0.19.0-rc.1.json) · [The five final MATLAB text changes and their source fingerprints](release-0.19.0-differences.json) · [Interval commands and result locations](../docs/interval_article.md). The public evidence summary omits absolute workstation paths. The original user-supplied report ZIPs are retained separately from the source distribution.

**Exact arithmetic has a separate record.** The native numerical runs did not request `Exact=true`; their exact-replay status is `NOT_REQUESTED`. They therefore do not establish native Symbolic execution of the rational certificates. The complete independent Python/SymPy replay was executed during development, as recorded below. The complete native MATLAB Symbolic route remains unverified.

## Auxiliary development evidence before the native run

The following checks were performed on 20 November 2026 before the native reports arrived. They remain separate evidence: auxiliary Octave runs and Python proof replay are not relabelled as native MATLAB runs.

| Scope | Actual evidence |
|---|---|
| Whole auxiliary core suite | GNU Octave 10.3.0: **646 passed, 0 failed, 85 skipped**, status `CORE_OK_NATIVE_PENDING`. The skipped cases require native MATLAB quaternion or the separately disabled graphics suite. |
| Complete interval numerical study | **24/24 groups, 110/110 comparisons** passed using ordinary component structs. The current-run archive includes separately recorded calculations and reference values. |
| New graphics | Three dedicated cases passed separately, checking eight figures and 32 PNG/PDF/FIG/CSV files, preserved caller state and unchanged saved numeric data. PNGs and rendered PDFs were inspected. This is auxiliary graphics evidence, not native MATLAB graphics acceptance. |
| Full independent exact replay | Python 3.12 / SymPy 1.14: 28 exact algebraic checks, 56 quartic Bernstein inequalities, all 24 quintic faces, 264 polynomial checks and 370 proof leaves at maximum depth five. |
| MATLAB-language symbolic replay | Octave symbolic 3.2.2 recomputed all 56 quartic inequalities. A separate quintic run verified the OOO centre and two complete OOO faces (22 polynomial checks), then was deliberately stopped because this bridge is slow. The complete quintic `.m` route and native Symbolic execution are not claimed tested. |
| Source integrity | The distribution manifest is regenerated from the final bytes with Python SHA-256. This Octave build has no Java backend, so its embedded MATLAB-style hash inventory explicitly reports unavailable hashes. Native `run_qstab_public` still requires its normal source and resolution checks. |

The auxiliary run initially exposed an existing Octave 10.3 parser issue when `Diagonal` also partially matched `Diagonalizable`. Registering the exact shorter name first fixes only option dispatch in `qstab_eigen_assess`; all numerical eigenpair checks are unchanged. All 74 affected legacy checks and the subsequent full suite pass. The baseline distribution itself was not edited.

[Machine-readable interval development evidence](interval-development.json) · [Full exact replay evidence](interval-exact-python.json) · [Reproduction commands](../docs/interval_article.md).

To reproduce native acceptance, run `run_qstab_public('Figures',false)` and `run_qstab_intervals('all','Plots',true)` in the normal public MATLAB environment. The optional `Exact=true` requests Symbolic Math Toolbox and records its result separately. Native tests use the actual MathWorks quaternion class.

## Historical validation of qstab 0.18.0

**Native MATLAB public acceptance: PASSED.** The release candidate `0.18.0-rc.1` was run on 9 November 2026 with MATLAB 24.2.0.2773142 (R2024b Update 2), Windows PCWIN64. The 693 individual regression checks all passed; none failed or were skipped. The Command Window transcript and the two archived outputs (`report.txt` and `report.mat`) identify the same run. The regression phase took 27.763 s; the complete public acceptance run took 39.130 s. These are measurements of this run, not benchmarks.

| What was checked | Evidence |
|---|---|
| Source identity | All 193 MATLAB source fingerprints and sizes recorded before and after the run match the 0.18.0-rc.1 distribution. MATLAB's search path and the caller's working folder were restored. |
| Clean public environment | The native MathWorks quaternion implementation was resolved. The static pre-run dependency scan of 40 public functions reported 99 files and no external files. This describes the branches and products identified by MATLAB; it is not a universal dependency proof. An external eigenvalue solver was not requested. |
| Native types, figures and guide examples | All 60 native quaternion/MATLAB cases and all 20 graphics cases passed. The guide test ran all eleven public keys. The run used `Figures=false`: ordinary graphics regressions (including export fixtures) ran; the optional publication-figure export/review was not requested. |
| Earlier evidence | Public 0.16.0: 661/661 on native R2024b Update 2. The candidate had already passed 613 core checks under GNU Octave 8.4.0 with 80 native/graphics cases skipped. [Historical Octave record](octave-core.json). These are separate runs. |
| Paper checks | Separately verified using the supplied exact arithmetic scripts: 47 mechanical checks; 15 four-generator subfamilies, 1,850 frequency cells and 18,500 Bernstein identities. These are not included in the count of 693. [Scope and replay instructions](../docs/papers.md). |

**From the tested candidate to 0.18.0.** Of 193 MATLAB files, 189 are byte-identical. Four carry updated version or descriptive text: `qstab_setup.m` (session version), `run_qstab_public.m` (version in an early-failure report), `src/+qstab/about.m` (version, date and validation text) and `src/+qstab/Contents.m` (help heading). No solver, tolerance, test, guide example or algorithm changed. The native report is attributed to its actual candidate version; a second native run of the renamed distribution is not claimed.

[Machine-readable native evidence](native-matlab.json) · [Candidate and final source fingerprints for the four changed files](release-differences.json). These records omit absolute workstation paths. The original native report ZIP was checked but is not part of the public distribution.

## Reproduce on another installation

From the extracted, writable package directory in a fresh MATLAB session:

```matlab
report = run_qstab_public('Figures',false);
disp(report.status)
disp(report.archiveFile)
assert(report.ok && report.publicAcceptance.ok, 'Read the public acceptance report.')
```

The archive, `qstab_test_results.zip`, contains report.mat and report.txt. Use `Figures=true` only to additionally export the publication figures. Report files from your own computer may contain installation paths. Numerical tests and plots exercise the implementation; they do not establish robust stability of an arbitrary new model.
