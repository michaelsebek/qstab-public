# Scalar one-parameter article: reproducibility data

The inputs correspond to **qB01 polynomial article, revision 7, 26 November
2026**. The manuscript is not redistributed in this folder. Its hashes and
the hashes of the supplied exact-verification source files are recorded in
`source_provenance.json`.

## MATLAB examples

After the usual public setup, run:

```matlab
R = run_qstab_oneparameter('all','Plots',true);
R.status
R.articleComplete
R.archiveFile
```

The numerical path uses MATLAB and its native `quaternion` object. With
`'UseQuaternion',false`, the same numerical algorithms receive ordinary
component structs instead. The optional exact path uses Symbolic Math
Toolbox; `ExactCertificates='auto'` checks availability and records whether
proofs were actually executed. To require them explicitly, use
`'ExactCertificates',true`. To run only the numerical comparisons, use
`'ExactCertificates',false`.

`R.ok` refers to the checks that ran. `R.articleComplete` additionally requires
all example groups and all exact proof groups. Numerical floating-point
bounds and unresolved tangencies remain labelled as such. The archive of a
run is always `<package>/qstab_oneparameter_results.zip`; the raw MAT, text,
CSV and figure files are under a unique subdirectory of
`<package>/oneparameter_results`.

`radius_results.csv` contains the numerical radius brackets, their status,
reference values, computed frequency values and exact rational bounds.
`ray_results.csv` distinguishes the signed numerical interval from the exact
enclosure of the **first positive** contact, and records the unnormalized
quaternion direction and original-polynomial lifting status. The general
check-by-check table is `comparison.csv`. A new ZIP is staged, unpacked and
compared byte for byte before replacing the previous result archive.

A complete exact run requires nine global-radius certificates, fourteen
first-ray proofs (including the added E5 variant), thirteen finite-contact
proofs on the original polynomial, and the 42 feedback/variant identities.
The fourteenth ray is unbounded and therefore has no contact to lift.

`inputs.json` contains the seven frequency datasets, their nine case/side
radius references, the thirteen original rational directional fixtures,
and the additional noncommuting E5 variant. Rational strings are preserved
separately from their numerical views. Dataset E7 uses the original large
coefficients; the MATLAB numerical runner records the explicit frequency
and amplitude normalization it applies.

## Independent exact Python replay

This is an explicit optional command, not a MATLAB dependency or a process
started by the toolbox:

```text
python -m pip install sympy==1.14.0
python reproducibility/qb01/replay_exact.py --output qb01_exact_results
```

It recomputes the original 13 directional cases (376 exact checks), 9
global radius certificates, and 42 additional identities for the mechanical
feedback interpretation and E5 variant. The scripts run in a temporary copy
and write a new report and logs in the chosen output directory. The packaged
historical outputs in `source/` are not overwritten. Their presence does not
count as a newly executed MATLAB or Python proof.

The exact radius certifier uses numerical reference values only to propose
rational endpoints and a rational frequency witness. Acceptance is decided
by fresh exact polynomial signs and root counts, not by the reference values.

## Provenance

The original files retained under `source/` are unmodified copies from the
author's supplied revision-7 reproducibility archive. Revision 7 itself was
editorial: the historical numerical benchmark and exact source files were
unchanged from revision 6. The standalone MATLAB runner and certificate
entry points are new qstab integrations. Recorded Python timings and SVD
witness errors are diagnostic observations, not universal test tolerances.
