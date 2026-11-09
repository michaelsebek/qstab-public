function info = about()
%ABOUT Version, scope and dependency information for qstab.
%
%   Examples (ascending coefficients):
%   I = qstab.about();
%   disp(I.dependencies);
%
%   Interpretation:
%   Version, numerical scope, dependency declarations and historical validation text; this is not a fresh validation run.
%
%   Reference: docs/reference/about.html.
%   See also qstab.about.
    info=struct();
    info.name='qstab - Quaternion Stability Toolbox';
    info.version='0.18.0';
    info.buildDate='2026-11-09';
    info.targetMATLAB='R2024b Update 2';
    info.inputTypes='MATLAB quaternion, numeric arrays, cells, and ordinary component/family structs';
    info.coefficients='Ascending: P(:,:,k+1) is the coefficient of v^k';
    info.nominalScope='Scalar algebraic zeros; right matrix spectra; Hurwitz/Schur/delay nominal tests';
    info.robustScope='Scalar zero exclusion for polytopes and correlated CT/DT families; full leading-hull bounds; qB02 cube/pair-free/quadratic and qB03 four-generator reductions';
    info.qB01='CT and discrete z/zinv scalar first-exit rays and one-quaternion global radius bounds; discrete state-matrix rays';
    info.coefficientBallScope='Full quaternion weighted l_eta coefficient balls; pointwise spectral classes and witnesses; CT/Schur/delay global radius brackets; no box/ball substitution';
    info.componentBoxScope='Actual-point intrinsic-rank zonotope gauges; CT/Schur/delay global brackets from nested complete-family tests; no sampled-direction lower-bound claim';
    info.experimental='Nominal regular square polynomial-matrix companion support only';
    info.dependencies={'Base MATLAB';...
        'A MathWorks product providing quaternion, when native quaternion inputs are used';...
        'Control System Toolbox optional for large nominal matrix Lyapunov solves'};
    info.rightEigenpairs='qstab.eigpairs: right eigenvectors, canonical diagonal values, per-column residuals and independent basis-rank diagnostics';
    info.validation='Public 0.18.0-rc.1: 693 passed, 0 failed, 0 skipped in native MATLAB R2024b Update 2. Final 0.18.0 changes release metadata and documentation only. See verification/VALIDATION.md.';
    info.benchmarkScope='Optional matched full/reduced scalar study; inputs, outputs and timing samples in report.benchmarks, unchanged two-entry ZIP';
    info.diagnosticScope='Read-only retained-cell coverage, fixed-argument outer sets, frequency-enclosure geometry and independently checked member witnesses';
    info.eigenComparison='Fourteen deterministic analytic spectra; optional explicit legacy handle, oriented eigenpairs and matched one-output timings; no qeig alias';
    info.paperExtensions='16 additional finite source-identity/member checks and eight selected regenerated illustrations; not full article reproduction';
    info.paperFigureDirectory='qstab_paper_figures in the package root; PDF/PNG/FIG/CSV outputs, not a replacement test-results ZIP';
    info.paperReproduction='20 source-linked selected qB01/qB02/qB03 cases with separate references and computed outputs; not all paper figures/examples';
    info.publicAcceptance='run_qstab_public: temporary package/MathWorks-only path, source and dependency audits, no legacy solver, unchanged two-entry result ZIP';
    info.releasePlan='Public release: tagged ordinary family structs, repeatable random generation and sampled zeros/values; no private class or solver required.';
    info.resultArchive='qstab_test_results.zip in the package root (beside test_qstab.m), independent of pwd and OutputDir';
    info.assurance='Floating point with residual/conditioning checks; no validated interval arithmetic';
    info.externalQuaternionPackagesRequired=false;
    info.customClasses=false;
    info.sourceProvenance='Installed .m SHA-256 inventory is embedded in the MAT report; explicit unavailable-hash fallback';
    info.wrapperPolicy='No private q-aliases or object adapters in public distribution. Separate integration package only.';
    if nargout==0
        fprintf('%s\nVersion %s (%s)\n',info.name,info.version,info.buildDate);
        fprintf('%s\n%s\n',info.nominalScope,info.robustScope);
        fprintf('%s\n',info.assurance);
    end
end
