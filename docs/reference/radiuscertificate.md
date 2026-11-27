# qstab.radiuscertificate

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf)

## Help and calling sequences

```text
RADIUSCERTIFICATE Verify rational lower/upper bounds for a quaternion radius.
  [bounds,Info] = qstab.radiuscertificate(F,rBounds)
  [bounds,Info] = qstab.radiuscertificate(P0,P1,rBounds)
  [bounds,Info] = qstab.radiuscertificate(P0,P1,rBounds,'Side','right')

  F is an unrestricted qstab.family('oneparameter',...) definition. P0 is
  a Hurwitz scalar polynomial and degree(P1)<degree(P0). One constant
  quaternion t forms P0+t*P1 for Side='left', or P0+P1*t for Side='right'.
  Coefficient vectors are ascending. Native quaternion values, numeric
  coefficients and ordinary component structs a,b,c,d are accepted.

  rBounds is [lower upper], or a two-element cell array of exact rational
  text, for example {'378937381963/1000000000000', ...
  '378937381964/1000000000000'}. The lower endpoint may be zero. Otherwise
  both endpoints must be finite, rational and 0<=lower<upper.

  bounds is the verified 1-by-2 enclosure, or [NaN NaN] if a requested
  proof fails. Inf Inf is returned for P1=0. Double endpoints are rounded
  outwards; Info.exactBounds contains the authoritative rational strings.
  Info.certified and Info.ok are true only when the nominal assumptions
  and BOTH bounds were proved in this invocation. Failed certificates
  are not instability decisions. No numerical search is silently used.

  The proof covers ALL frequencies, using x=omega^2 and exact rational
  Sturm root counts of Q_r=delta0-2*r^2*delta01+r^4*delta1. A positive lower
  endpoint is accepted if Q_lower has no nonnegative real root. For the
  upper endpoint, a root of Q_upper proves that a critical gain of norm
  no greater than upper exists. Alternatively UpperX supplies one exact
  nonnegative rational x at which Q_upper<=0 or T_upper=n0-upper^2*n1<=0.
  Repeated roots and flat extrema are included by square-free reduction.

  Options:
    Side            'left' (default) or 'right'.
    RationalData    [] (default), or struct('P0',C0,'P1',C1), where C0,C1
                    are 4-by-N ascending real rational component matrices.
                    Entries may be numeric, symbolic rational, or fraction
                    strings in cells. Explicit data must round to the
                    supplied input coefficients, without shape changes.
    UpperX          []: count roots of Q_upper. Otherwise certify the
                    supplied x=omega^2 directly; it is NOT omega.
    Variable        's' only. Grade must be [].
    ExtendParameter false. A finite family domain requires explicit true;
                    the result then concerns its law on all quaternion t.

  Exact numeric semantics: a supplied double denotes its exact binary
  value (sym(value,'f')); no rational approximation is guessed. Use
  fraction text in RationalData when the intended source is, for example,
  exactly 1/10 rather than its stored binary approximation. Symbolic input
  coefficients must be rational. Returned reports contain no sym objects.

  Example (requires Symbolic Math Toolbox):
    [b,I] = qstab.radiuscertificate([1 1 1],1,{'0.8660','0.8661'});
    assert(I.certified); b

  Requirement: optional Symbolic Math Toolbox. Ordinary numerical calls
  to qstab.stabilityradius/freqmargin remain independent of that toolbox.
  This function does not modify files, paths, symbolic assumptions or RNG.
  See also qstab.raycertificate, qstab.stabilityradius, qstab.freqmargin.
```
