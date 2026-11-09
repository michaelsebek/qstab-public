# qstab · Quaternion Stability Toolbox for MATLAB

**Version 0.18.0** · [Read the illustrated guide](docs/qstab_manual.pdf) · [Function reference](docs/api_index.md) · [Installation and dependencies](docs/installation.md)

Test nominal and robust strict stability for scalar quaternionic polynomials; compute right spectra of quaternion state matrices and representatives of algebraic zero classes. Families may be polytopes, independent coefficient intervals, real affine families, or one constant quaternion acting from a specified side. Coefficient balls and component boxes, continuous and discrete tests, radii, witnesses, and sampled plots are included in this release. Robust polynomial matrices are outside the supported scope; the square polynomial-matrix companion is a nominal experimental feature.

**Native MATLAB verification:** all 693 regressions passed on R2024b Update 2 (Windows) using the release candidate. The final package changes version information and documentation only. [What was verified](verification/VALIDATION.md).

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

The initial companion manuscript is M. Šebek, *A Sharp Face Theorem for Robust Stability of Quaternionic Polynomial Polytopes*, prepared for Automatica. The public package also includes methods not covered by that paper. [Research notes and reproduction key](docs/papers.md) · [Cite this version](CITATION.cff) · [MIT license](LICENSE) · [Release notes](RELEASE_NOTES.md).

For optional integrations and MathWorks product checks, see [installation](docs/installation.md). There is no mandatory third-party GitHub download.
