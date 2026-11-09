% qstab - Quaternion Stability Toolbox, public release 0.18.0.
%
% Family definitions and exploration (standard MATLAB structs)
%   family           - Define interval, polytope, affine or one-quaternion uncertainty.
%   randfamily       - Draw a random family of a chosen kind, degree and variable.
%   plotFamilyZeros  - Sample algebraic zero classes in the representative plane.
%   plotFamilyValues - Sample quaternion values in the class/slice or norm plane.
%
% Public stability and spectrum functions
%   isstable         - Root-free nominal tests, with an explicit spectral option.
%   robuststable     - Whole-family Lyapunov/zero-exclusion methods.
%   zeroexclude      - qB02/qB03 norm-value exclusion; optional four-generator reduction.
%   roots            - Algebraic zero class representatives, including infinity.
%   eig              - Canonical right eigenvalue class representatives.
%   raymargin        - CT and discrete first exit in a specified direction.
%   freqmargin       - Scalar critical distance at CT frequency or DT angle.
%   stabilityradius  - Scalar global CT or unit-circle parameter-radius bounds.
%
% Visualization
%   plotSpectrum     - Exact 2-D quotient and 3-D slice of similarity classes.
%   valueset         - Reduced norm-polynomial values, with correlated samples.
%   plotValueSet     - Display value samples and zero, without a certificate claim.
%
%   valueenclosure   - Separate fixed-argument outer polygons/rectangles.
%   plotEnclosure    - Outer sets and samples without a false certificate claim.
%   certificateData - Read retained accepted/unresolved cell records.
%   certificateEnclosure - Reconstruct one stored whole-frequency enclosure.
%   plotCertificate - Frequency projections or correctly labelled slack.
%   witnessData     - Original-family membership and independent spectral check.
%   plotWitness     - Input/forward algebraic classes of a checked member.
%
% Numerical data and operators
%   components       - Native quaternion/numeric/component-struct parts.
%   polycomponents   - Scalar ascending coefficient components and grade.
%   familycomponents - Real affine descriptions of scalar uncertainty families.
%   chi              - Complex adjoint representation.
%   normpoly         - Real central norm polynomial, ascending coefficients.
%   companion        - Normalized matrix or generalized companion pencil.
%   lyapoperator     - Real Lyapunov operator on quaternion Hermitian matrices.
%   about            - Version, dependencies, and scope.

% qB02/qB03 extensions
%   leadingbound    - Whole-family leading-coefficient separation.
%   kharitonovcube  - Exact CT eight-parameter interval reduction.
%   zeroexclude    - Joint affine/discrete intervals and interval-face reductions.

% Reproducible study (examples/):
%   run_qstab_benchmarks - Matched full/reduced comparisons, return-only.
%   test_qstab('Benchmarks',true) saves the study in the usual report ZIP.

% Diagnostic study (examples/):
%   run_qstab_diagnostics - Data-only output and optional separate figures.
%   test_qstab('DiagnosticExamples',true) retains data in the unchanged ZIP.

% Weighted scalar coefficient balls (qB02):
%   ballvalue  - Pointwise minimum radius to an algebraic zero class.
%   ballset    - Guarded membership at supplied spectral query points.
%   ballradius - Global CT/Schur/delay radius bracket, including degree loss.
%   ballstable - Strict robust stability at a prescribed ball radius.

% Discrete single-parameter examples (examples/):
%   run_qb01_discrete - Schur/delay rays, shared-parameter radii, degree loss.
%   test_qstab('DiscreteExamples',true) retains these data in the same ZIP.

% Componentwise scalar boxes (qB02):
%   boxvalue   - Actual quaternion point gauge, primal/dual bounds and witness.
%   boxstable  - Complete component-box test at a specified scale.
%   boxradius  - Global scale brackets from whole-family tests and witnesses.
%   run_qb02_boxes - Component-box examples, including a hidden edge margin.
%   test_qstab('BoxExamples',true) retains data in the unchanged report ZIP.

% Source-linked selected examples (in the examples directory):
%   qstab_paper_cases  - catalogue of 20 source-linked fixed cases.
%   run_qstab_papers   - selected-source checks; no writes or plotting.

% Additional article support (examples folder):
%   run_qstab_paper_extensions - Sixteen additional finite source checks.
%   qstab_paper_figures        - Eight ordinary-data mathematical illustrations.
%   plot_qstab_paper_figure    - Render one saved item, without solver calls.
%   export_qstab_paper_figures - PDF/PNG/FIG/CSV export from saved data.

% Right-spectrum study (examples/):
%   qstab_eigen_cases          - Fourteen deterministic analytic reference matrices.
%   run_qstab_eigen_comparison - Optional explicit existing qeig comparison, no source edits.
%   test_qstab('EigenComparison',true) saves ordinary data in the same report ZIP.

% Right eigenvectors (separate from the established eig diagnostics API):
%   eigpairs - [V,D,info], with per-column residuals and quaternionic rank.

% Standalone acceptance (package root, no new numerical methods):
%   qstab_check_installation - Read-only source integrity and name resolution.
%   run_qstab_public         - Isolated PUBLIC acceptance; same result ZIP.
