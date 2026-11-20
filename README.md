# qstab · Quaternion Stability Toolbox for MATLAB

**Version 0.19.0** · [Read the illustrated guide](docs/qstab_manual.pdf) · [Function reference](docs/api_index.md) · [Interval article: reproduce all examples](docs/interval_article.md) · [Installation and dependencies](docs/installation.md)

Test nominal and robust strict stability for scalar quaternionic polynomials; compute right spectra of quaternion state matrices and representatives of algebraic zero classes. Families may be polytopes, independent coefficient intervals, real affine families, or one constant quaternion acting from a specified side. Coefficient balls and component boxes, continuous and discrete tests, radii, witnesses, and sampled plots are included in this release. Robust polynomial matrices are outside the supported scope; the square polynomial-matrix companion is a nominal experimental feature.

**Validation:** the 0.19.0-rc.1 sources passed **731 native MATLAB regression checks**, with no failures or skips, on R2024b Update 2. The separate interval study matched **all 24 example groups and 110 comparisons** and exported eight figures. The final 0.19.0 changes only version/documentation text and two plot labels; numerical routines and tests retain their validated contents. Optional exact-rational replay has a separate scope and evidence record. [What was verified](verification/VALIDATION.md).

### Reproduce the interval article and supplement

```matlab
qstab_setup('Examples',true);
R = run_qstab_intervals('all','Plots',true);
R.archiveFile
assert(R.ok,'Inspect the saved interval report: %s',R.archiveFile)
```

This runs the 24 source-linked groups from `IJRCN_manuscript_v7` and `IJRCN_Zenodo-supplement_v5`. Use `run_qstab_intervals('list')` to see the key or `run_qstab_intervals('I19','Plots',true)` for the five-gain regulator. Each run retains actual results, separate reference values, comparison checks and diagnostics. The current-run ZIP is `qstab_interval_results.zip` beside the launcher. Native quaternion is the default; `UseQuaternion=false` explicitly selects ordinary component structs. The optional exact certificate replay is separate from numerical root checks; see [instructions and scope](docs/interval_article.md).

### First minute

Unzip to a new folder. Change the MATLAB current folder to the directory containing `qstab_setup.m`. Run:

```matlab
qstab_setup('Examples',true);
F = qstab.family('polytope',{[1 1],[2 1]});    % p(s)=s+a, 1 <= a <= 2
[tf,info] = qstab.robuststable(F);
tf                                             % 1
info.status
R = qstab_guide_example('G06');                % one checked, headless example
```

Use `tf==1`, `tf==0`, and `isnan(tf)`: the three results are strictly stable, demonstrably not strictly stable, and undecided or unmet assumptions. Open `info.status` for the distinction. For a robust failure also inspect its witness. Coefficients are **ascending**; [2 3 1] is 2+3s+s². Scalar roots are algebraic similarity classes: they need not vanish under point evaluation at the displayed representative.

### From question to runnable code

| Question | Where to start |
|---|---|
| How do I enter quaternion coefficients or select continuous, forward, or delay time? | Guide examples G01 and G02 |
| What are the right eigenvalues of a matrix? | G03 · `qstab.eig`, `qstab.eigpairs` |
| How far can one quaternion move in a fixed direction? What if its physical range is bounded? | G04 and G11 · `qstab.raymargin` |
| How do independent component intervals differ from shared parameters? | G05 · `qstab.family` |
| How can a local enclosure exclude zero when the global one cannot? | G06 and G07 · [polytope article and exact checks](docs/papers.md) |
| How do I inspect a certificate or a destabilising member? | G08 |
| What changes when the perturbation set is a ball instead of a box? | G09 |
| How do I plot sampled zero classes or the values of a seeded random family? | G10 · `qstab.plotFamilyZeros`, `qstab.plotFamilyValues` |

`qstab_guide_example()` lists the keys; `qstab_guide_example('all')` executes all eleven in local function workspaces. The examples do not write files, download software or open figures. To make figures, use the explicit plot calls in the guide. To validate a full **native** installation separately, use `report=run_qstab_public('Figures',false)`. That acceptance run writes `qstab_test_results.zip` beside the runner; it checks basic graphics even when optional figure export is disabled. [Read the validation record](verification/VALIDATION.md).

### Article, citation and license

The companion texts include M. Šebek's preprint [*A Sharp Face Theorem for Robust Stability of Quaternionic Polynomial Polytopes*](https://doi.org/10.5281/zenodo.24136448), submitted to Automatica, and the manuscript *Kharitonov-Type Face Reductions and Exact Stability Tests for Quaternionic Interval Polynomials*. The published technical note [*Quaternionic Interval Polynomials: Worked Examples and Verification Guide*](https://doi.org/10.5281/zenodo.24389775) supplies worked examples and verification guidance. The software version and the mathematical texts have separate citation records. The package also includes methods outside those texts' scope. [Polytope research notes](docs/papers.md) · [Interval study and references](docs/interval_article.md) · [Cite this version](CITATION.cff) · [MIT license](LICENSE) · [Release notes](RELEASE_NOTES.md).

For optional integrations and MathWorks product checks, see [installation](docs/installation.md). There is no mandatory third-party GitHub download.
