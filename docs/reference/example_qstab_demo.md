# qstab_demo

[Function index](../api_index.md) · [Manual](../qstab_manual.pdf) · [MATLAB source](../../examples/qstab_demo.m)

## Help and calling sequences

```text
QSTAB_DEMO Small examples of quaternion nominal and robust stability.
  OUT = qstab_demo uses the native MATLAB quaternion class and opens
  spectral/value-set figures. The scalar examples use ASCENDING
  coefficient tensors P(:,:,k+1)=P_k; matrices are ordinary quaternion
  arrays. No custom quaternion or polynomial class is needed.

  OUT = qstab_demo('Plots',false) runs without graphics.
  OUT = qstab_demo('UseQuaternion',false,'Plots',false) uses ordinary
  a,b,c,d component structs. This path can be used for an auxiliary
  Octave check; it does not validate MATLAB's native quaternion class.

  Demonstrations:
    p(s)=s^2+3*s+2+j: nominal root-free Hurwitz test and algebraic zeros;
    an upper-triangular quaternion matrix: Lyapunov test and right spectrum;
    z-(0.4+0.2*j) and 1-(0.4+0.2*j)*z^-1: reciprocal conventions;
    scalar first-degree quaternion intervals and a three-vertex polytope;
    optional exact 2-D/3-D similarity-class views and sampled value sets.

  The robustness examples have simple independent analytic explanations:
  every member is s+a0 with Re(a0)>=1, hence its zero class is Hurwitz.
  Numerical robustness outputs still carry their ordinary floating-point
  assurance; the plotted value samples are not a robustness certificate.
```
