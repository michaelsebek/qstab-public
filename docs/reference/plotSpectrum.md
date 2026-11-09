# qstab.plotSpectrum

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../src/+qstab/plotSpectrum.m)

## Help and calling sequences

```text
PLOTSPECTRUM Plot right-spectral or algebraic quaternion similarity classes.
  [AX,Info] = qstab.plotSpectrum(lambda, ...)
  AX is the target MATLAB axes. Info contains all previous plot handles
  and the plotted data or diagnostics. Request Info for reuse or inspection.
  Info is the full graphics-handle structure returned by earlier releases.
  [~,H] = qstab.plotSpectrum(LAMBDA) displays each class by its canonical
  complex representative a+1i*rho, rho = norm(Im(q)) >= 0. Numeric complex,
  MATLAB quaternion, and component-struct input accepted by
  qstab.components are supported. Multiplicities are retained in H.

  Options:
    'Variable'   's' (default), 'z', or 'zinv'. Stability regions are
                 Re(lambda)<0, abs(lambda)<1, and abs(lambda)>1.
    'Dimension'  2 (default) or 3.
    'Parent'     Existing axes. Their contents and hold state survive.
    'Visible'    'on' (default) or 'off' for a newly created figure.
    'Tolerance'  Nonnegative relative boundary tolerance; default 100*eps.

  The 2-D plot is the exact quotient of quaternion similarity classes.
  The 3-D plot is the EXACT SLICE Im_k(q)=0, not an arbitrary projection:

      (Re(q), Im_i(q), Im_j(q)) = (a, rho*cos(phi), rho*sin(phi)).

  Thus a nonreal class becomes a circle, and a real class becomes a point.
  The continuous boundary is the plane Re(q)=0. The discrete boundary is
  the unit sphere. An entire similarity class has one stability status.
  For 'zinv', the stable side is OUTSIDE the sphere.

  Infinite values are counted and annotated separately, never moved to
  finite plotting coordinates. Their system interpretation requires the
  polynomial grade/descriptor convention and is not inferred here.
  The plot classifies supplied values; it is not a robustness certificate
  and does not establish completeness of a computed spectrum.

  Examples:
      qstab.plotSpectrum([-1+2i; .3+.7i],'Variable','s');
      qstab.plotSpectrum([-1+2i; .3+.7i],'Variable','z','Dimension',3);
  See also qstab.eig, qstab.roots, qstab.components, qstab.plotValueSet.

  Examples (ascending coefficients):
  [lambda,I] = qstab.roots([2 3 1]);
  [ax,H] = qstab.plotSpectrum(lambda,'Variable','s','Dimension',3);

  Interpretation:
  2D uses the quotient (real part, imaginary norm). 3D uses a norm- and real-part-preserving slice of a similarity class, not an arbitrary projection of H. Infinite representatives are identified separately.

  See also qstab.about.
```

## Additional option notes

| Option | Default | Meaning |
|---|---|---|
| `variable` | `'s'` | Time convention: s (continuous), z (forward shift) or zinv (delay). Reinterprets the supplied coefficients; it does not transform a system. |
| `dimension` | `2` | 2 for canonical complex representatives, 3 for a real-part/norm-preserving class slice. |
| `parent` | `[]` | Existing axes handle; its current content and hold/NextPlot state are preserved by the public plotting routine. |
| `visible` | `'on'` | on/off for a newly created figure. Does not alter the data or the validity of a result. |
| `tolerance` | `100*eps` | Numerical residual/separation or pairing guard as described for this routine. It is not a physical uncertainty radius and not a directed-rounding proof. |
