# Install qstab and check what MATLAB will run

Extract the complete release into a new, writable folder. The package root contains `qstab_setup.m`, `src/+qstab` and `examples`. GitHub's “Download ZIP” may give the outer folder a different name; select the directory containing `qstab_setup.m`. Do not merge release trees, add `+qstab` directly, or use `addpath(genpath(...))`: old code and same-named functions can then shadow the intended ones.

In MATLAB select that folder and run:

```matlab
qstab_setup('Examples',true);
which qstab.isstable -all
S = qstab_guide_example('G01'); assert(S.ok);
```

Setup adds the root, `src` and optionally `examples` for this session. It does not save the path, edit startup, install toolbox products, or add global aliases. The guide uses ordinary numeric arrays or structs with equally sized real `a,b,c,d` fields, so all its numbered examples can be tried without the MathWorks quaternion class.

| Task | What to install or connect | How to check |
|---|---|---|
| Numeric and component-array public APIs; examples G01–G11 | Base MATLAB; no external GitHub functions. | `qstab_guide_example('all')`. |
| Native quaternion input/output and the full native acceptance run | A licensed MathWorks product providing [quaternion](https://www.mathworks.com/help/fusion/ref/quaternion.html), e.g., Sensor Fusion and Tracking Toolbox. One provider is enough. Use [MathWorks Add-Ons](https://www.mathworks.com/help/matlab/matlab_env/get-add-ons.html) or the MathWorks installer for your release, then restart MATLAB if requested. | `which quaternion -all`; the selected constructor should be under `matlabroot`. Try `q=quaternion(1,2,3,4); [a,b,c,d]=parts(q);` |
| Large nominal matrix certificate requiring a Lyapunov/Stein solve beyond the configured Kronecker limit | Optional Control System Toolbox, providing [lyap](https://www.mathworks.com/help/control/ref/lyap.html) and [dlyap](https://www.mathworks.com/help/control/ref/dlyap.html). Install through MathWorks if you need this branch. Smaller solves are included. | `which lyap -all`, `which dlyap -all`; inspect the diagnostic if the optional solve is unavailable. |
| Compare an existing *right* eigensolver against qstab | Optional; only for a user-requested comparison. Install the other solver following its own README, put only its entry directory on the session path, inspect its name with `which ... -all`, and pass an explicit trusted function handle. | For an installed compatible `qeig`: `E=run_qstab_eigen_comparison('LegacyFunction',@qeig,'UseQuaternion',true,'Timing',false);`. Do this separately from public acceptance. It must accept a native matrix and return right eigenvalues, and optionally [V,D]. |
| Replay the manuscript's exact arithmetic | Python 3.8+ standard library; no Python packages or MATLAB toolboxes. | [Run the two article checks](papers.md). |

[qpe-public](https://github.com/michaelsebek/qpe-public) is a separate toolbox for **evaluation zeros**. If you also need that mathematical problem, follow its own installation instructions in another folder. It is not an implementation of `qstab.roots`. Neither a left-eigenvalue solver nor a third-party quaternion matrix multiplication routine is required by qstab.

**Full check:** Start a fresh MATLAB session, select the qstab folder, and use `report=run_qstab_public('Figures',false); assert(report.ok && report.publicAcceptance.ok);`. This runner temporarily restricts resolution to qstab and MATLAB installation paths, checks the source manifest and static dependencies and restores the path and current directory on exit. Use a writable package root: it writes fresh report files and refreshes the fixed `qstab_test_results.zip`. Back up an earlier report if you wish to keep it. Temporary path isolation is not a sandbox for untrusted code. No add-ons are installed or downloaded by qstab. See [test evidence](../verification/VALIDATION.md).

[Home](../README.md) · [Manual](qstab_manual.pdf) · [All entry points](api_index.md).
