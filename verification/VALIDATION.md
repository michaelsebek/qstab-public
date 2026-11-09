# Validation of qstab 0.18.0

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
