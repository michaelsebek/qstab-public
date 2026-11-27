# qstab.raycertificate

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf)

## Help and calling sequences

```text
RAYCERTIFICATE Prove the first positive loss on one quaternion ray exactly.
  [gammaBounds,Info] = qstab.raycertificate(F,d)
  [gammaBounds,Info] = qstab.raycertificate(P0,P1,d)
  [gammaBounds,Info] = qstab.raycertificate(P0,P1,d,'Side','right')

  F is an unrestricted qstab.family('oneparameter',...) definition. The
  raw form describes P0+t*P1 (Side='left', default) or P0+P1*t ('right'),
  with t=gamma*d and gamma>0. d is one NONZERO quaternion and is never
  normalized. P0 must be Hurwitz, degree(P0)>=1 and degree(P1)<degree(P0).
  Coefficients are ascending numeric/native quaternion vectors or ordinary
  structs of the four real components a,b,c,d. Only continuous time is used.

  gammaBounds=[lo hi] encloses the FIRST POSITIVE real gamma where strict
  stability is lost, including isolated tangential contacts. It is not
  the signed [negative positive] interval returned by qstab.raymargin.
  Call with -d for the opposite direction. [Inf Inf] proves no positive
  contact. A finite root outside the isolation budget gives [NaN NaN].
  Exact rational endpoints are Info.exactBounds; returned doubles are
  rounded outwards. The gain norm is gamma*abs(d), not gamma itself.

  The real norm polynomial is quadratic in gamma. Its reduced Hurwitz
  determinant is formed exactly, then its square-free part is isolated
  with exact Sturm counts. Even-multiplicity roots are retained. No
  floating-point eigenvalue is accepted as a proof or rounded to real.
  The nominal Hurwitz property and fixed degree are also checked exactly.

  Options:
    Side            'left' (default) or 'right'.
    RationalData    [] or struct('P0',C0,'P1',C1,'Direction',d4).
                    C0,C1 are 4-by-N ascending rational component arrays;
                    d4 is 4-by-1. Numeric, sym rational and cell fraction
                    strings are accepted. Direction may be omitted.
                    Explicit data must round to the supplied coefficients.
    AbsTol          1e-12, requested absolute gamma enclosure width.
    MaxBisections   256, positive integer. Exact counting is independent
                    of this precision budget; inspect widthTargetMet.
    VerifyLifting   true: verify the contact on the original polynomial
                    where an exact algebraic representation is available.
                    Lifting status is separate from the exact first-root
                    proof and is never replaced by a numerical residual.
    Variable        's' only. Grade must be [].
    ExtendParameter false; finite family bounds require explicit true.

  Info.certified/ok is true only after an exact first-root proof (or an
  exact proof of no positive root). Info.isolation stores the determinant,
  square-free polynomial, Sturm sequence and root counts. Reports contain
  ordinary MATLAB data and rational strings, so they can be saved normally.
  Numeric doubles are interpreted as exact stored binary rationals using
  sym(value,'f'); decimal/fraction source data belong in RationalData.

  Example (a tangential contact, requires Symbolic Math Toolbox):
    P0=struct('a',[1 2 1],'b',[0 0 0],'c',[1 0 0],'d',[0 0 0]);
    P1=struct('a',[-1 -1],'b',[0 0],'c',[-1 -1],'d',[0 0]);
    [ab,I]=qstab.raycertificate(P0,P1,1); assert(I.certified); ab
    % The exact first positive contact is gamma=1, despite re-entry.

  Optional dependency: Symbolic Math Toolbox. Numerical raymargin remains
  available without it. Symbolic determinant computations may be costly.
  See also qstab.raymargin, qstab.radiuscertificate, qstab.freqmargin.
```
