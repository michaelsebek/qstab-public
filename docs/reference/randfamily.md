# qstab.randfamily

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/randfamily.m)

## Help and calling sequences

```text
RANDFAMILY Draw a repeatable, well-defined scalar quaternion polynomial family.
  F = qstab.randfamily                         random interval, degree 1, 2 or 3
  F = qstab.randfamily(D)                      interval of degree D
  F = qstab.randfamily(D,KIND)                 choose one of the four family kinds
  F = qstab.randfamily(D,KIND,N)               N polynomials or affine directions
  F = qstab.randfamily(KIND)                   choose the kind; draw the degree
  F = qstab.randfamily(KIND,D,...)             equivalent kind-first spelling
  F = qstab.randfamily(...,'int',...,'Seed',42) small integer coefficients, repeatable

  D is the degree in the forward variable s or z. [] chooses uniformly
  from 1:3. In delay notation the coefficients are reversed at grade D.
  D=0 generates nonzero quaternion constants, except for 'oneparameter',
  which requires a nonzero perturbation of strictly lower forward degree.

  KIND           Family and meaning of N
  'interval'     Default. Independent bounds on the real,i,j,k parts of
                 every free coefficient. No N.
  'oneparameter' One quaternion t shared by all powers. P0+t*P1 (left)
                 or P0+P1*t (right). No N. The nominal is built to be stable.
  'polytope'     Real convex hull of N defining polynomials (default 3).
                 N counts generators, not necessarily extreme points.
  'affine'       P0+sum(theta_j*G_j) with N shared real parameters
                 (default 2). One interval per parameter.

  Optional inputs
  's','z','d' or 'Variable',V  Time variable. Aliases: p for s, q for z,
       zi/z^-1/zinv for delay. Default: s.
  'Width',W    Positive coefficient scale (default 1). Scales the interval
       half-widths, the displacements from the first convex generator, or
       the perturbation directions. It is not a stability radius.
  'int'        Integer defining coefficients; W must be a positive integer.
       This does not discretize the parameters, weights or family members.
  'Side',S     'left' (default) or 'right'; only for 'oneparameter'.
  'Bounds',B   'oneparameter': finite ordered 4-by-2 bounds on real,i,j,k
       of t, default [-1,1] for each. [] means all of H; use for margins.
       'affine': finite ordered N-by-2 bounds; each defaults to [-1,1].
       Equal endpoints fix that coordinate. Bounds are independent of W.
  'SamplingPeriod',H   0 in continuous time; positive in discrete time
       (discrete default 1, in the same units as your model).
  'Seed',K     Integer 0..2^32-1. Repeats the family without advancing
       the caller's random stream. Without K the active stream advances.
  'Degree',D, 'Kind',KIND (or 'Type',KIND) and 'Count',N are also accepted.
       Use each argument once. N is a positive integer for polytope/affine.

  F is an ordinary definition struct accepted by qstab.robuststable,
  qstab.raymargin (oneparameter), and qstab.plotFamilyValues / plotFamilyZeros.
  Its defining polynomials are component structs with row vectors a,b,c,d.
  For D>=1 every member is monic in the forward variable. Delay coefficients
  are reversed at the full grade; the current-sample coefficient is one.
  The nominal is stable by construction in exact arithmetic. No stability
  test of a whole random family is implied.

  F=qstab.randfamily(3,'interval','Seed',17);
  [~,V]=qstab.plotFamilyValues(F,1,9,'Plot',false);
  F=qstab.randfamily(2,'polytope',4,'z','Seed',2);
  [tf,I]=qstab.robuststable(F,'MaxCells',16);
  F=qstab.randfamily(2,'oneparameter','Bounds',[],'Seed',17);
  ab=qstab.raymargin(F,1);
  See also qstab.family, qstab.familycomponents, qstab.plotFamilyZeros.
```
