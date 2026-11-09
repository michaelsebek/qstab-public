# qstab.leadingbound

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/leadingbound.m)

## Help and calling sequences

```text
LEADINGBOUND Nonvanishing of the FULL leading-coefficient uncertainty set.
  A qstab.family definition supplies its variable, grade and parameter model.
  [LOWERBOUND,Info] = qstab.leadingbound(data,...)
  LOWERBOUND is the computed lower bound on the required coefficient norm.
  A zero bound alone does not prove degree loss; inspect Info.status.
  Info is the full diagnostic structure returned by earlier releases.
[~,R] = qstab.leadingbound(V,'Variable','s') accepts the same family formats as
familycomponents. For zinv the current-sample coefficient is tested.
R.status is 'separated', 'zeroWitness', or 'undetermined'; R.lowerBound is
a sufficient lower bound on the quaternion norm over the entire family.
A small minimum-norm residual is NOT an assumption-failure certificate.
Supporting directions are found by a bounded Frank-Wolfe iteration, but
acceptance uses a directly computed support minimum. No optimizer toolbox
or quaternion class is needed for component inputs. validated is false.

  Examples (ascending coefficients):
  [lowerBound,R] = qstab.leadingbound({[1 1],[2 2]},'Variable','s');
  disp(lowerBound);

  Interpretation:
  A nonzero leading coefficient must be excluded over the whole original family, not just its generators. For delay notation the relevant coefficient is checked after grade-preserving reversal.

  See also qstab.about.
```
