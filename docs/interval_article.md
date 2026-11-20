# Reproduce the interval article and supplement

This chapter accompanies **qstab 0.19.0**, prepared from the numerically validated 0.19.0-rc.1 candidate. This final packaging is not a fresh native MATLAB acceptance run. It covers every concretely specified numerical example in *Kharitonov-Type Face Reductions and Exact Stability Tests for Quaternionic Interval Polynomials*, manuscript v7, and *Quaternionic Interval Polynomials: Worked Examples and Verification Guide*, supplement v5. General existence and complexity theorems without numerical inputs are not replaced by invented fixtures.

The public companion supplement is [*Quaternionic Interval Polynomials: Worked Examples and Verification Guide*, version 1.0](https://doi.org/10.5281/zenodo.24389775), published as a technical note on 20 November 2026. The precise manuscript v7 and supplement-source v5 revisions used by this study are identified by hashes below; source revision labels and public record versions are distinct. The main interval manuscript has no journal DOI assigned here.

## Start with one command

Run these commands in the directory containing `qstab_setup.m`, in the ordinary PUBLIC MATLAB environment:

```matlab
qstab_setup('Examples',true);
report = run_qstab_intervals('all','Plots',true);
assert(report.ok,'Inspect %s',report.archiveFile);
report.archiveFile
```

The study computes all 24 groups, keeps the calculated values and paper references separately, and creates a new `qstab_interval_results.zip` **in the package root**. It creates that archive even if a numerical example fails. Detailed results are in `report.mat`; readable comparisons are in `report.txt` and `comparison.csv`. Each run has a unique raw-results subdirectory, given by `report.outputDirectory`. Only the latest run enters the fixed-name ZIP, so old figures cannot silently enter a new report.

Native MathWorks `quaternion` is the default input format. Base MATLAB can execute the same calculations with ordinary four-component structures:

```matlab
report = run_qstab_intervals('UseQuaternion',false,'Plots',true);
```

This second route is explicitly labelled as component input. It does not establish native quaternion integration. No private classes or external quaternion package are required.

To list or select examples:

```matlab
catalogue = run_qstab_intervals('list');
R = run_qstab_intervals('I19');
R.cases(1).output.radius

Q = run_qstab_intervals({'I09','I10'},'Plots',false);
Q.cases(1).output.quadraticForm
```

For an in-memory run with no files or plots, use `'Save',false`. Example keys remain stable when printed article numbering changes. A selected run stores its results in `cases(1)`, `cases(2)`, and so on, in catalogue order. For the complete run, I19 is `report.cases(19)`.

## Figures and saved data

`Plots=true` exports eight figures for a complete run, each as PNG, PDF, FIG and CSV: the running square, its three norm-value segments, the Hurwitz and Schur boundary counterexamples, the edge reserve, the five-gain gap and three state trajectories, and the two-gain frequency reserve. Every plot states whether it displays an analytical region or sampled curve. Existing figures and axes are preserved.

To replot a saved calculation without rerunning the solvers:

```matlab
saved = load(fullfile(report.outputDirectory,'report.mat'),'report');
files = plot_qstab_intervals(saved.report, ...
    'Directory',fullfile(report.outputDirectory,'replotted'),'Visible','off');
```

Only computed, selected cases are plotted. The external replot command does not rewrite the original result ZIP. The full `run_qstab_intervals(...,'Plots',true)` call saves its plots in the ZIP automatically.

## What is being compared?

A quaternion polynomial is entered with ascending coefficients and four real component rows. The real norm is

\[
N_p(s)=p_0(s)^2+p_1(s)^2+p_2(s)^2+p_3(s)^2.
\]

The norm's complex zeros encode algebraic classes. A complex representative need not itself be an actual quaternion evaluation zero. The study records both representations, and explicitly evaluates the claimed quaternion right roots in the original polynomial. It preserves tiny nonzero coefficients such as the quintics' `epsilon=1e-12`.

Every comparison has a name, a calculated `actual` value, an independent `reference`, a tolerance, and a `pass` flag. `report.ok` states whether those selected comparisons passed. It does not turn a finite root plot into a continuum certificate. The numerical status returned by `qstab.isstable` is retained even if a separate exact calculation settles the example. In particular, default tolerances may return `NaN` for extremely scaled quintics; the replay does not weaken those guards globally.

## Coverage map

| Key | Calculation | Article / supplement |
|---|---|---|
| I01 | Algebraic representative versus actual right zero | Article §2 |
| I02 | Correlation between real norm coefficients | Article §2 |
| I03 | Running square, exact value ranges, corner failure, both edge contacts, correlated line | Article running example; supplement §3.1-3.4 |
| I04 | Declared degree, zero at infinity, loss of regularity | Article §2 |
| I05 | Stable endpoints with an unstable one-component interior | Article one-component theorem |
| I06 | Hurwitz square with stable complete boundary and unstable centre | Article quadratic sharpness; supplement §4.1 |
| I07 | One-component slice and reserve `3-sqrt(5)` | Supplement §3.5 |
| I08 | Actual-point gauges and unreachable directions | Supplement §3.6 |
| I09 | OOO quartic, simple-root derivatives, boundary quadratic form | Article §5.1; supplement §6.1 |
| I10 | EOO quartic and its different contact scale | Article §5.2; supplement §6.2 |
| I11-I14 | Four explicit degree-five parity classes | Article exact degree-five witnesses |
| I15 | Isolated angular minimum at a rank drop, reserve `40/11` | Article global reserves; supplement §7.2 |
| I16 | Nonmonic Schur cube, degrees 1, 2, 5 and extra zero modes | Article Schur sharpness |
| I17 | Monic Schur quadratic and annular stability region | Article monic quadratic; supplement §4.2 |
| I18 | Monic Schur cubic and radial Schur inequalities | Article monic cubic |
| I19 | Five independent regulator gains, exact reserve equation, 3×401 trajectory samples | Article regulation application |
| I20 | Two independent gains in one component, 402 frequency samples | Article restricted application |
| I21 | Edge-interior contact missed by every vertex ray | Supplement §7.3 |
| I22 | Zero-frequency, degree-loss and infinite reserves | Supplement §7.4 |
| I23 | Full linear and quadratic coefficient boxes | Supplement §8.3 |
| I24 | Degree/parity continuous-face and algebraic-candidate counts | Article combinatorial tables |

The source anchors are also stored with each result. The script does not modify either manuscript.

## Three new reusable functions

### Independent coefficient uncertainty in one component

`qstab.boxfreqmargin` considers independent perturbations in one Hamilton component of the polynomial coefficients. For example, changing real damping and real stiffness independently is a two-real-parameter uncertainty, not one shared quaternion parameter.

```matlab
P = quaternion([1 1 1],[2 0 0],[0 2 0],[0 2 0]);
W = zeros(4,3);
W(1,1:2) = 0.1;
[r,Info] = qstab.boxfreqmargin(P,W,0.7,'Witness',true);
```

At frequency omega the uncertain component ranges over a complex rectangle. The other three components determine two complex exclusion points. The required radius is the smaller of the two rectangular gauges:

\[
\mu(\omega)=\min_{\sigma\in\{-1,1\}}
\max\left\{\frac{|\operatorname{Re}(\sigma\zeta-\widehat p_\nu)|}{W^e},
\frac{|\operatorname{Im}(\sigma\zeta-\widehat p_\nu)|}{W^o}\right\},
\qquad
\zeta^2=-\sum_{\kappa\ne\nu}\widehat p_\kappa^2.
\]

Zero parity width fixes that coordinate exactly. The implementation does not project a small nonzero mismatch to zero. `Info` records both branches and, if requested, an actual original-coefficient polynomial and quaternion right-root witness. This is a result at the given frequency. A general global minimum over frequency remains a separate task. The general multi-component actual-point gauge is still `qstab.boxvalue`.

### Continuous faces appropriate to the degree and parity

```matlab
F = qstab.family('interval',L,U);
[faces,Info] = qstab.intervalfaces(F);
[tf,detail] = qstab.zeroexclude(F,'Reduction','interval-faces', ...
                              'MaxSubfamilies',1536);
```

`faces` contains complete affine families, not selected nominal members. The combined test enforces degree preservation before a stability decision. For full reduced cubes the counts are 256 in degree one, 1536 in degrees two and three, 896 in degree four, and 1024 from degree five. Lower-dimensional cases with inactive directions are handled explicitly. The legacy `kharitonov-faces` selection remains unchanged; the specialized monic `quadratic-faces` method can still be cheaper, with at most 192 families. These counts describe continuous families, not individual evaluations or root solves.

### Local motion of an algebraic zero

```matlab
[gradient,Hessian,Info] = qstab.rootmotion(F,1i,'Parameters',theta);
Q = -real(Hessian)/2;
```

The function computes derivatives of a **simple** zero of the norm polynomial by implicit differentiation. This reproduces both quartic contact matrices without finite differences. It returns unresolved data for an unverified or multiple norm zero. A Hessian describes local geometry; the finite-box proofs below additionally control the remainder and the complete boundary.

A complete small example is \(p(s,t)=s+1+t+\mathbf j\), with \(-0.1\le t\le0.1\). The upper representative is \(\lambda(t)=-1-t+i\), so its gradient is -1 and its Hessian is zero:

```matlab
P = quaternion([1 1],[0 0],[1 0],[0 0]);
F = qstab.family('affine',P,{1},[-0.1 0.1]);
[g,H,Info] = qstab.rootmotion(F,-1+1i);
g                         % -1
H                         % 0
Info.valid                % true
```

MATLAB's `1i` denotes the ordinary complex representative here; the polynomial's native quaternion constant contains the Hamilton unit j.

## Recompute the exact continuum certificates

The default study uses floating-point MATLAB arithmetic. The article's exact quartic and quintic certificates are also included as optional reproducible calculations.

With Symbolic Math Toolbox:

```matlab
proof = verify_qstab_interval_certificates('Part','all');
assert(proof.ok)
proof.coverage
```

Alternatively, include this fresh proof in the normal result archive:

```matlab
report = run_qstab_intervals('all','Exact',true,'Plots',true);
assert(report.ok,'Inspect %s',report.archiveFile)
```

The exact function reconstructs rational coefficients; it does not convert rounded double coefficients to symbolic values. It recalculates all 56 quartic tensor-Bernstein inequalities, the four central quintic Routh columns, and the complete closed boundaries of all four quintics. Their 24 faces require 264 polynomial positivity checks, covered by 370 exact dyadic proof rectangles. Coverage is checked, so a missing cell cannot count as a complete proof. The tenth Hurwitz determinant is the positive constant coefficient times the ninth.

The independent Python/SymPy route in `reproducibility/qb02` additionally checks exact scalar identities and reserve isolators. It has been executed independently during development. It is optional and is not called by the MATLAB numerical functions. The MATLAB symbolic route still needs its first native acceptance run.

These are concrete certificates for the supplied examples. They do not implement the paper's entire general algebraic elimination algorithm or general algebraic global-radius solver. `boxradius` retains its existing global brackets and can stop with a gap. A dense angular mesh cannot repair the isolated rank-drop example I15; its exceptional directions must be included analytically.

## Source identity and publication status

The published starting point audited for this work was GitHub tag `v0.18.0`, commit `21a8d39fca6dff14273e6e2cb825ae108db45501`, archived at [doi:10.5281/zenodo.24086085](https://doi.org/10.5281/zenodo.24086085). That DOI identifies the older software version. The present 0.19.0 package is prepared locally; these files do not claim that a new GitHub release or Zenodo software version has already been published. Paper and supplement identifiers are separate from a software-version identifier.

Exact manuscript source SHA-256:

```text
IJRCN_manuscript_v7.tex
0c95bdc5b8c948b13ca38e4cde762ce74f1b529eaeb96a0fce9d2a8916e18098
IJRCN_Zenodo-supplement_v5.tex
664b23dd602b3aca58d2826f49afd857e7899f9d628e03e6da60ba044a197a19
```

On 20 November 2026, the **0.19.0-rc.1 numerical baseline** passed 731/731 regression checks in MATLAB R2024b Update 2. The separate native interval run matched 24/24 groups and 110/110 comparisons and exported all eight figure groups. Its exact replay was `NOT_REQUESTED`. The 0.19.0 final documentation and two label-position adjustments were prepared afterward; no new native run of that final packaging is claimed. The MATLAB symbolic path still needs its first native validation, and private adapters have their own separate validation. See [the validation ledger](../verification/VALIDATION.md) for the original run identities and scope.
