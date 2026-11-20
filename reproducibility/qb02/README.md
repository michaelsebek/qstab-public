# Exact checks for the interval-polynomial article

These optional independent verifiers reproduce the main article's explicit
quartic certificates, the four degree-five boundary counterexamples, their
central Routh counts, and selected exact identities and reserve isolators.
They complement the ordinary numerical MATLAB study.

## Run the exact calculation

The independent Python route needs **Python 3.10 or later and SymPy 1.12 or
later**. It does not need MATLAB, NumPy, SciPy, a numerical root solver, or an
internet connection after SymPy is installed.

```text
python -m pip install sympy
python reproducibility/qb02/replay.py --output my_exact_interval_report.json
```

Run the second command from the public package root, or use an absolute path
to `replay.py`. The output report is a freshly computed result; it identifies
its actual Python and SymPy versions and explicitly sets `nativeMATLABRun`
to false. Without `--output`, no report file is created.

An independent MATLAB route is provided for installations with **Symbolic
Math Toolbox**:

```matlab
proof = verify_qstab_interval_certificates('Part','all');
assert(proof.ok)
proof.coverage
```

Use `'Part','quartic'` or `'Part','quintic'` to run one part. The MATLAB route
reconstructs the exact original polynomial coefficients and recomputes the
Bernstein bounds and, for quintics, the Hurwitz determinants. It does not call
Python and does not use stored pass flags as proof. It returns a report and
creates no files. Its first native MATLAB execution remains to be checked by
the maintainer; the Python replay has already been independently executed.
Symbolic determinant calculation may take longer than the numerical study.
Neither optional dependency is required for the ordinary public functions.

## What is proved

- **28 exact algebraic checks:** quartic norm factorizations and cubic norm
  Routh columns, four quintic central Routh columns and their rational bounds,
  OOO root-bracket/shift inequalities, the running example and one-component
  re-entry determinants, the monic Schur cubic inequalities, and exact
  reserve isolators and identities.
- **56 strict rational tensor-Bernstein inequalities** reproduce the OOO and
  EOO quartic continuum estimates. Every rational Bernstein coefficient is
  recalculated from the stated polynomials, including the full shift interval.
- **24 complete closed quintic faces** cover all six boundary faces of each
  OOO, EOO, EEE and EEO parameter cube. On each two-dimensional face, the
  degree-ten real norm is constructed from the original quaternion components.
  Its leading and constant coefficients and its first nine Hurwitz
  determinants are strictly positive on the entire square. The tenth
  determinant equals its constant coefficient times the ninth.
- These **264 polynomial positivity checks** use **370 exact dyadic proof
  leaves**, at maximum subdivision depth five. Every leaf contains a whole
  rectangle. Coverage and exact positivity are replayed, so neither the
  continuous faces nor the shift intervals are replaced by sampled points.

The centres' exact Routh sign sequences contain two sign changes. Consequently
the four norm polynomials each have exactly two open right-half-plane roots.
The whole stable boundary and the unstable centre are separate conclusions.
The coefficient `epsilon = 1/10^12` is preserved exactly in all calculations.

## Files and trust boundary

- `verify_exact_examples.py` derives the exact algebraic and quartic data;
  `exact_example_certificate.json` retains their rational results.
- `verify_quintic_boundaries.py` derives the norm polynomials and Hurwitz
  determinants; `quintic_boundary_certificate.json` retains coefficient hashes
  and a complete dyadic partition with exact rational lower bounds.
- `replay.py` checks SHA-256 hashes in `sources.json`, then executes both
  independent verifiers. Recalculated mathematical evidence, not a stored
  `ok` flag, determines success.

To regenerate the proof partitions, run each verifier with `--build`. This
rewrites its certificate and is intended for development; a changed
certificate also needs a deliberately updated source manifest before the
combined replay will accept it.

The reference text is `IJRCN_manuscript_v7.tex`, supplied for this software
audit. Its SHA-256 is recorded in `sources.json`; the unpublished manuscript
itself is not redistributed in the software. The certificates provide a
finite sufficient proof for these explicit examples. They are **not an
implementation of the article's general exact fixed-box decision procedure
or general global algebraic reserve algorithm**. Ordinary qstab numerical
answers retain their documented floating-point assurance level.
