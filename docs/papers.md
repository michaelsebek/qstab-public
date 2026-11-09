# Companion paper and reproducible code

**Michael Šebek, “A Sharp Face Theorem for Robust Stability of Quaternionic Polynomial Polytopes.”** Manuscript prepared for *Automatica*, source version 6 November 2026. No acceptance, journal DOI or preprint identifier is claimed. Cite the software version separately using [CITATION.cff](../CITATION.cff).

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

**Other methods.** Interval, one-quaternion, coefficient-ball and matrix routines are documented and runnable in the same release; the polytope paper does not establish all their claims. See the [manual](qstab_manual.pdf) for assumptions, budgets and what each return value means.
