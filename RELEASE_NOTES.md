# qstab 0.20.0

Adds the scalar one-quaternion study, exact first-positive-ray certificates (`qstab.raycertificate`) and verification of supplied global-radius bounds (`qstab.radiuscertificate`). Extends fixed-frequency gain geometry, ordered SVD diagnostics, original-polynomial evaluation witnesses and frequency scaling. Decimal/scientific text and polynomial zero identities now remain exact throughout the symbolic branch. See [the one-parameter guide](docs/oneparameter_article.md).

Native MATLAB R2024b Update 2 evidence is attributed to the actual versions: the rc.1 transcript reports 748/748 ordinary public checks; rc.2 passed 32/32 targeted checks and 13/13 article groups with 139/139 comparisons, including 9 radius certificates, 14 ray certificates (13 finite contacts and one unbounded ray), 13 finite original-polynomial liftings and 42 exact identities. The 221 installed MATLAB source hashes in the article report match the delivered rc.2. Numerical tangencies can remain unresolved while exact first-contact proofs succeed; the report retains both outcomes.

Relative to the tested rc.2, finalization changes four MATLAB files only in version/descriptive strings, plus documentation and publication metadata. The other 217 MATLAB files, all numerical and symbolic algorithms, all examples and all test implementations are byte-identical. The complete 748-check suite was not rerun on rc.2; no new native execution of the renamed final distribution is claimed. [Validation](verification/VALIDATION.md) records the earlier failed exact rc.1 run and the successful rc.2 replacement separately.

Exact proofs require optional Symbolic Math Toolbox and apply to continuous-time scalar polynomial families with rational component data and fixed leading coefficient. The global-radius API verifies user-supplied bounds. It does not claim a general exact radius algorithm for uncertain quaternion matrices.

---

# qstab 0.19.0 · 20 November 2026

Adds a reproducible study of the interval manuscript v7 and supplement v5, with 24 source-linked example groups, independent comparison values, eight figure exports, and an automatically rebuilt result ZIP. New public functions are `qstab.intervalfaces`, `qstab.boxfreqmargin`, and `qstab.rootmotion`. The opt-in `Reduction='interval-faces'` uses the new continuous degree/parity-aware selection in zero exclusion; existing reduction modes retain their previous meaning.

The optional exact replay reconstructs rational polynomials and checks complete quartic inequalities and quintic boundary faces. It is not a sample-based proof and does not claim implementation of the paper's general exact fixed-box or global-reserve algorithms. MATLAB symbolic replay is an optional layer; standard numerical calls acquire no Symbolic dependency. Tiny quintic coefficients are retained, and unresolved floating-point stability flags are not overridden by the example reference values.

**Native validation passed:** the 0.19.0-rc.1 sources ran on MATLAB R2024b Update 2, Windows PCWIN64, on 20 November 2026. All **731 regression checks passed**, with zero failures and zero skips: 646 core, 62 native-type and 23 graphics checks. All 208 installed MATLAB source fingerprints and byte sizes matched the delivered candidate before and after the run. A separate native interval-study run matched **24/24 example groups and 110/110 comparisons**, and saved all eight figures with their numerical plot data.

**From the tested candidate to 0.19.0:** four MATLAB files change only version or descriptive text (`qstab_setup.m`, `run_qstab_public.m`, `src/+qstab/about.m`, and `src/+qstab/Contents.m`). Two labels in `examples/plot_qstab_intervals.m` are corrected for readable display; the plotted values are unchanged. Numerical algorithms, tolerances, example calculations and regression tests are unchanged. Documentation and citation metadata now describe the completed validation. These records identify the actual tested candidate; they do not claim another native execution after the final text edits.

The optional native MATLAB Symbolic replay was not requested by the numerical acceptance runs. Its complete native execution remains unverified; the independent Python/SymPy exact replay has its own recorded evidence. See [the validation ledger](verification/VALIDATION.md) and [interval reproduction instructions](docs/interval_article.md).

## qstab 0.18.0 · 9 November 2026

Based on the complete public 0.17.0 source. This release makes the family constructor, random family generation, directional domain reporting and sampled plots available to standalone users. New family values are MATLAB structs. The core conventions stay the same: ascending coefficients, result first and diagnostics second, canonical representatives of algebraic zero classes and right spectra. The existing positional calls remain available. The MATLAB guide uses a fixed list of 11 executable keys.

The companion manuscript concerns polytopes. Other implemented methods are available in the same package; a software release does not assert that each method has been published. Future paper references can be added to CITATION.cff and [the reading notes](docs/papers.md). The optional article replay is separate from the MATLAB solvers.

Native public acceptance for 0.18.0-rc.1 passed on MATLAB R2024b Update 2 (Windows) on 9 November 2026: 693 passed, zero failed and none skipped. The reported digests of all 193 MATLAB files match the candidate archive. Between that archive and 0.18.0 only the version displays, the reporting text, a help heading and distribution documentation/metadata changed. All numerical routines, all examples and all tests retain their tested contents.

[Validation record](verification/VALIDATION.md) · [Installation](docs/installation.md) · [License](LICENSE).
