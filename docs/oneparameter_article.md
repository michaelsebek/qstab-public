# Scalar one-parameter article: calculations, certificates and saved results

Version 0.20.0 accompanies Michael Šebek's *Stability Radii of Quaternionic Polynomial Families With One Constant Quaternion Parameter*, supplied revision 7 of 26 November 2026. Submission to *Systems & Control Letters* is planned. This description does not assert journal acceptance or a newly assigned DOI. The manuscript itself is not redistributed.

## Run one example or the whole article

From the installed public root:

```matlab
qstab_setup('Examples',true);
run_qstab_oneparameter('list');
A = run_qstab_oneparameter('E1','Plots',true);
A.status
A.archiveFile

A = run_qstab_oneparameter('all','Plots',true, ...
    'ExactCertificates','auto');
A.status
A.articleComplete
A.archiveFile
assert(A.ok,'Inspect the article report: %s',A.archiveFile)
```

The first two calls use the same catalogue as the full run. `Cases={'E1','E5','E5V'}` selects several keys. Defaults use MATLAB's native `quaternion` inputs. `UseQuaternion=false` explicitly selects ordinary real component structures. MATLAB never starts Python or installs a product during these calls.

| Key | Source question | What is actually checked |
|---|---|---|
| E1 | Left versus right calibration | Both fixed-frequency formulas, ordered critical gains, nonunit rays and all-frequency radius bounds. |
| E2 | Singular perturbation frequency matrix | A rank-one perturbation, a finite contact, its original-polynomial witness and ray. |
| E3 | Multiple nearest gains | A sphere with center −1/2 and radius 1/√2; every shown critical parameter has norm √3/2. |
| E4 | Near-degenerate gain geometry | A unique or unresolved sensitive gain; a tiny singular-value gap is not a proof of a sphere. |
| E5 | Tangential first loss | The first contact at one, even if the numerical ray remains unresolved and stability later returns. |
| E5V | No common complex slice | The additional noncommuting tangency at one. |
| E5C | Norm-induced multiplicity | The real cubic's first contact at one, without skipping an even-multiplicity determinant zero. |
| E6 | A flat optimum | The small critical frequency and distinction between a numerical bracket and an exact certificate. |
| E7 | Physical scaling | E1 with amplitude 10^80 and frequency scale 10^6; parameter radii are unchanged. |
| D11–D13 | Directional diagnostics | No finite contact, rejection of a later nonimaginary guardian candidate, and spherical contact. |
| A01 | Assumptions and conventions | Zero perturbation, nonmonic coefficients and an unnormalized direction. |

The catalogue includes nine case-side global radius problems and fourteen rays, including the new nonslice tangency. It preserves the article's rational source data and reference values separately from actual function outputs. Machine-dependent runtimes and SVD witness errors remain diagnostics rather than portable pass constants.

## Outputs and evidence levels

`A.outputDirectory` is a unique run folder under `<public root>/oneparameter_results`. It contains the MAT report, readable text, comparison CSV and requested figures. `A.archiveFile` is always `<public root>/qstab_oneparameter_results.zip`. The ZIP is staged, unpacked and compared with the current files before replacement. Earlier raw run folders remain available; the fixed ZIP contains the most recent invocation, including a single-example invocation.

- `A.ok` means that the executed comparisons passed and explicit requirements were met.
- `A.articleComplete` additionally requires the complete catalogue, all required fresh exact certificates and the 42 additional identities.
- `A.numericalAllResolved` tracks the numerical ray results independently. A difficult tangency can correctly remain unresolved in the numerical branch while an exact first-root proof succeeds.
- `A.exact` states whether exact arithmetic was actually executed. An archived source certificate is never silently counted as a fresh proof.
- Input hashes, source provenance, resolved function paths and current source inventory remain in the report. The unsigned manifest detects consistency, not publisher identity.

Errors inside individual examples are recorded before the final archive is made. Invalid options or an unavailable explicitly requested native input facility fail at preflight. Use a writable package root. The runner restores its temporary MATLAB path additions and does not overwrite caller variables or initialized symbols.

## Numerical margins and gain geometry

For the scalar CT model

$$
p_L(s,t)=p_0(s)+t p_1(s),\qquad p_R(s,t)=p_0(s)+p_1(s)t,
$$

one quaternion is shared by every coefficient. The multiplication side belongs to the model. The nominal polynomial must be Hurwitz and the perturbation must have lower degree for these CT radius methods.

`qstab.freqmargin` minimizes the norm of the parameter over the imaginary quaternion sphere at a prescribed frequency. `qstab.stabilityradius` additionally covers all frequencies and returns a numerical lower/upper bracket. `qstab.raymargin` instead fixes a direction and returns the negative and positive endpoints of the connected stable interval around zero. With `t=alpha*u`, endpoint `alpha` has parameter norm `abs(alpha)*abs(u)`; `u` is not normalized.

The new `Info.geometry{k}` describes a nearest-gain point or sphere, or explicitly unresolved/empty geometry. For a sphere, its ordinary real description is

$$
t=c+r_s Bv,\quad B\in\mathbb R^{4\times3},\quad B^TB=I,\quad \|v\|_2=1.
$$

Rows are the real, i, j and k coordinates. `r_s` is the sphere's radius about its center, not the norm radius returned by the first output. `Info.zeroWitnesses{k}` records right- and left-evaluation zeros on the original critical polynomial, their status and residuals. A `sphericalCandidate` needs exact confirmation; a complex class representative need not be an evaluation zero at that orientation.

An independent ordered 2×2 complex SVD comparison is reported without a new Control System Toolbox dependency. It is a floating-point diagnostic, not a proof of exact singular-value multiplicity or a generic quaternion-matrix stability-radius formula.

CT scalar margins accept `FrequencyScale='auto'` or a positive scalar. They may work with `s=tau*sigma` and a common positive coefficient scaling. Parameter radii, multiplication order and ray coordinates are unchanged; frequencies and boundary points are returned in the original units. `Info.scaling` records the transformation. The nominal `isstable` decision rule is not weakened.

## Optional exact methods

These calls require **Symbolic Math Toolbox**. Ordinary numerical functions do not.

```matlab
[provedBox,I] = qstab.radiuscertificate([1 1 1],1, ...
    {'0.8660','0.8661'});
assert(I.certified)
provedBox
I.exactBounds
```

`radiuscertificate` verifies supplied rational bounds over the whole frequency half-line. It does not automatically optimize a new bracket. Both the lower and upper proof must succeed. Failure to verify proposed bounds is not an instability decision.

```matlab
P0 = struct('a',[1 2 1],'b',[0 0 0], ...
            'c',[1 0 0],'d',[0 0 0]);
P1 = struct('a',[-1 -1],'b',[0 0],'c',[-1 -1],'d',[0 0]);
F = qstab.family('oneparameter',P0,P1);
[gBox,I] = qstab.raycertificate(F,1);
gBox                 % first positive contact: [1 1]
I.certified
I.lifting
```

`raycertificate` isolates the first **positive** ray contact using an exact square-free Hurwitz determinant and Sturm counts. Multiple roots and tangencies are retained. This output is an enclosure of one nonnegative number, not the signed interval from `raymargin`. Call with the opposite direction for the negative side. `[Inf Inf]` proves no positive contact; unresolved isolation returns NaNs. `AbsTol` and `MaxBisections` control refinement; a proof can be valid even if its requested width was not reached, which is reported by `widthTargetMet`.

`Info.certified` records the first-contact proof. `Info.lifting.exact` separately records an exact verification on the original quaternion polynomial. The default `VerifyLifting=true` requests it; numerical residuals are not a replacement for that proof. Implicit algebraic descriptions may be retained rather than untrustworthy decimal roots.

Stored doubles denote their exact binary rationals. Decimal text such as `'0.1'`, scientific text such as `'1e-1'`, and fraction text such as `'1/10'` denote the same exact rational number. Candidate rc.2 constructs these values from integer digits and powers of ten; it does not rely on a symbolic floating-point text conversion. A symbolic input must itself be an exact integer/rational object; `vpa` approximations are rejected. To certify intended decimal polynomial coefficients such as precisely 1/10, provide `RationalData` with 4-by-N ascending component arrays of rational strings. They must round to the supplied coefficients and preserve exact zeros; no rational approximation is guessed. Exact endpoint strings are authoritative; double output endpoints are rounded outwards. Reports contain ordinary serializable MATLAB data, not symbolic objects.

For a bounded family, exact functions require `ExtendParameter=true` and state that their result concerns the unrestricted model. They are scalar continuous-time methods, with `Variable='s'` and no declared extra grade. General uncertain quaternion state matrices and discrete exact certificates are outside this new scope.

## Regression checks and independent replay

```matlab
R = run_qstab_public('Figures',false);
assert(R.ok && R.publicAcceptance.ok)
T = test_qstab_oneparameter('ExactCertificates','auto');
T.status
T.complete
T.archiveFile
```

The main acceptance checks native quaternion, graphics, source resolution and numerical regressions. The targeted suite additionally exercises the reusable exact methods when available. Their archives are `qstab_test_results.zip` and `qstab_oneparameter_test_results.zip`. Missing optional exact work is not labelled complete.

To require every article proof:

```matlab
A = run_qstab_oneparameter('all','Plots',true,'ExactCertificates',true);
assert(A.articleComplete,'Incomplete proofs: %s',A.archiveFile)
```

An independent Python/SymPy replay is available as an explicit command from the package root:

```text
python reproducibility/qb01/replay_exact.py
```

See [its README](../reproducibility/qb01/README.md) for dependencies and outputs. This replays the preserved original thirteen-ray/376-check calculations, nine all-frequency certificates and 42 additional identities. It is separate evidence from executing the new MATLAB APIs. The supplied article and replay source hashes are in [source_provenance.json](../reproducibility/qb01/source_provenance.json).

Native MATLAB R2024b Update 2 reports for rc.2 record 32/32 targeted checks and the complete 13-group article study with 139/139 comparisons. Fresh exact work includes 9 global-radius certificates, 14 ray certificates (13 finite contacts and one unbounded ray), 13 finite original-polynomial liftings and 42 identities. All 221 MATLAB source hashes recorded in the article report match the supplied rc.2. The ordinary numerical ray results for E5, E6, E5V and E5C remain unresolved; their separate exact certificates succeed. This is retained as an intended distinction in the report.

The final 0.20.0 preserves those numerical, symbolic, example and test implementations and updates version/descriptive text, documentation and metadata. The earlier 748/748 ordinary public result belongs to rc.1; its separate 27/28 exact result remains a failed historical run. No new full public or final-version native execution is claimed. See [the validation ledger](../verification/VALIDATION.md) and [public verification commands](TESTING_CZ.txt).
