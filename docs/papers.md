# Companion papers and reproducible code

## Polytope preprint

**Michael Šebek, “A Sharp Face Theorem for Robust Stability of Quaternionic Polynomial Polytopes.”** Preprint version 1.0, published on 11 November 2026 and submitted to *Automatica*: [doi:10.5281/zenodo.24136448](https://doi.org/10.5281/zenodo.24136448). This identifier belongs to the preprint. The included reproducibility material retains the supplied manuscript source version of 6 November 2026, as recorded in [the source inventory](../reproducibility/qb03/sources.json). Cite the software version separately using [CITATION.cff](../CITATION.cff).

A shared vector of real nonnegative weights defines every coefficient of a polytope member. Testing the generators alone is inadequate. In the scalar constant-degree setting covered by the paper the sharp face result reduces the tests to whole faces of dimension at most three; a four-generator implementation covers complete subfamilies, including interiors. The package also offers direct family tests. A successful numerical subdivision and a verification in exact rational arithmetic are different forms of evidence.

| Key | Reproduce | What it shows |
|---|---|---|
| G06 | `qstab_guide_example('G06')` | Two-generator teaching model: dividing the segment excludes zero at one frequency. |
| G07 | `qstab_guide_example('G07')` | Six-generator cubic of the manuscript at zero frequency: a zero of the outer hull, but positive values for the selected actual members. |
| Mechanical exact check | See command below | Quaternion identities, determinants, rational inequalities and the independent-box witness for the four-mode mechanical example. |
| Cubic exact replay | Same command | An independent checker reconstructs the stored rational certificate: 15 four-generator subfamilies, 1,850 frequency cells, 18,500 Bernstein identities. The local lower bound 13730490575/1073741824 is for the tested faces on [0,30]; the analytical all-frequency bound 1296/361 has a different scope. |

From the repository root, if you want the paper checks in addition to MATLAB:

```text
python reproducibility/qb03/replay.py
```

Use `python3` on systems where `python` does not select Python 3. The dispatcher verifies the seven article files against their source digests, copies them into a temporary directory, and runs the two supplied checkers there with assertions enabled. No installation, network request or mutation of the distributed certificate is needed. The included generator is available for research but is not invoked by the replay. The source files reside in [reproducibility/qb03](../reproducibility/qb03/). The manuscript itself is not redistributed here.

## Interval manuscript and technical note

The interval study uses Michael Šebek's manuscript **“Kharitonov-Type Face Reductions and Exact Stability Tests for Quaternionic Interval Polynomials”**, source version 7, together with the supplied supplement version 5. The main article remains a manuscript with its own publication record.

The companion technical note **“Quaternionic Interval Polynomials: Worked Examples and Verification Guide”**, version 1.0, was published on 20 November 2026: [doi:10.5281/zenodo.24389775](https://doi.org/10.5281/zenodo.24389775). This DOI identifies the technical note. The precise source fingerprints used by the software study are listed in the [interval reproduction chapter](interval_article.md#source-identity-and-publication-status).

Run `run_qstab_intervals('all','Plots',true)` after `qstab_setup('Examples',true)` to compute all 24 source-linked groups, retain their 110 result comparisons and export eight figures with their numerical data. See [the interval reproduction chapter](interval_article.md) for individual example keys, output locations, dependencies and the separate optional exact-rational replay. Ordinary MATLAB numerical decisions and exact certificate replay provide distinct forms of evidence; their executed scopes are recorded in [the validation ledger](../verification/VALIDATION.md).

**Other methods.** One-quaternion, coefficient-ball and matrix routines are documented and runnable in the same release. Each method retains its own assumptions and references. See the [manual](qstab_manual.pdf) for assumptions, budgets and what each return value means.
