# qstab 0.18.0 · 9 November 2026

Based on the complete public 0.17.0 source. This release makes the family constructor, random family generation, directional domain reporting and sampled plots available to standalone users. New family values are MATLAB structs. The core conventions stay the same: ascending coefficients, result first and diagnostics second, canonical representatives of algebraic zero classes and right spectra. The existing positional calls remain available. The MATLAB guide uses a fixed list of 11 executable keys.

The companion manuscript concerns polytopes. Other implemented methods are available in the same package; a software release does not assert that each method has been published. Future paper references can be added to CITATION.cff and [the reading notes](docs/papers.md). The optional article replay is separate from the MATLAB solvers.

Native public acceptance for 0.18.0-rc.1 passed on MATLAB R2024b Update 2 (Windows) on 9 November 2026: 693 passed, zero failed and none skipped. The reported digests of all 193 MATLAB files match the candidate archive. Between that archive and 0.18.0 only the version displays, the reporting text, a help heading and distribution documentation/metadata changed. All numerical routines, all examples and all tests retain their tested contents.

[Validation record](verification/VALIDATION.md) · [Installation](docs/installation.md) · [License](LICENSE).
