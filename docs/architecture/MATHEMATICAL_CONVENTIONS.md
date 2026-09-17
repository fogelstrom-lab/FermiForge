# Mathematical conventions for the spinful 3He solver

Status: draft 0.1 for implementation review. Statements labelled **source**
transcribe the paper or current code. Statements labelled **project** are
proposed conventions for the new 2D interfaces. An unresolved item must not be
silently converted into a project convention.

## 1. Source hierarchy and names

The primary physics notation is `dcvlong.pdf`, especially Eqs. (1)-(3),
(7)-(13), (21)-(22), and Appendix A. The behavioral implementation references
are `incoming/legacy-f77` and `new_src`. When paper and code differ, an adapter
must name the conversion and a regression test must cover it.

The paper distinguishes the physical order parameter `A` from the
constant-spin-rotation-removed field `A_tilde` through

```text
A_{alpha i}(r,phi) = Delta_0 sum_beta R_{alpha beta} A_tilde_{beta i}(r,phi).
```

The code arrays are named `dxx` through `dzz`; they do not encode in their type
whether the constant rotation and bulk amplitude have been removed. Until this
is verified against the generating calculation, documentation must call them
the **code-basis A field**, not silently identify them with either paper field.

## 2. Coordinates and indices

**Project:** use right-handed Cartesian coordinates `(x,y,z)`, with `z` along
the vortex/cylinder axis and translation invariance `partial_z=0` in the first
2D implementation. Cylindrical coordinates obey
`x=r cos(phi)`, `y=r sin(phi)`, and
`e_phi=(-sin(phi),cos(phi),0)`.

**Source:** `A_{alpha i}` is a complex 3 by 3 matrix. The first index
`alpha in {x,y,z}` is spin and the second index `i in {x,y,z}` is orbital.
The code storage and `op_xyz` output order are

```text
Axx Axy Axz  Ayx Ayy Ayz  Azx Azy Azz,
```

with real and imaginary columns interleaved after the radius. **Project:** the
2D field API shall store this logical order explicitly as
`A(spin,orbital,point)`; memory-layout choices belong to the execution layer.

Spatial rotations act on the orbital index, spin rotations on the spin index,
and a gauge transformation multiplies every component by the same `exp(i chi)`.
These actions must remain separate APIs.

## 3. Triplet gap matrix

**Source (paper Eq. 1):** for momentum direction `qhat=p/pF`,

```text
Delta_hat(r,qhat) = sum_{alpha,i} A_{alpha i}(r)
                    i sigma_alpha sigma_y qhat_i.
```

Equivalently, with `d_alpha=sum_i A_{alpha i} qhat_i`,
`Delta_hat=(d dot sigma) i sigma_y`. This gives

```text
Delta_hat = [ -d_x + i d_y,   d_z
                d_z,          d_x + i d_y ].
```

The `new_src/riccati.f90` construction of its 2 by 2 coherence matrix has this
same component pattern when the scalar/singlet component `gamma(0)` vanishes.
This exact formula, rather than a remembered sign convention, is the proposed
boundary between `A` storage and spin-Nambu algebra.

Pauli matrices have their usual ordering `(sigma_x,sigma_y,sigma_z)` and the
spin basis is `(up,down)`. Any superconducting singlet term is to be added as
`Delta_s i sigma_y`; it is zero in the present 3He calculation.

## 4. Gauge, circulation, and axial harmonics

**Source:** a singly quantized far field satisfies
`A_tilde_{alpha i} -> delta_{alpha i} exp(i phi)`. The circulation number is
`M=1`; the input variable `vort` supplies this winding in the code through
`exp(i*vort*phi)`.

The complex cylindrical basis uses spin projection `sigma=-1,0,+1` and orbital
projection `ell=-1,0,+1`. Paper Eq. (8) assigns an angular factor

```text
exp(i Q phi) exp(i (M-sigma-ell) phi).
```

Only `Q=0` occurs in the axial solver; even nonzero `Q` is allowed by the
double-core `2m'm'` symmetry. The code's harmonic output order is

```text
A++ A+0 A+-  A0+ A00 A0-  A-+ A-0 A--.
```

The Cartesian-to-harmonic formulas currently implemented in `init_calc.f90`
are authoritative for regression. For example,

```text
A++ = (Axx-Ayy-i(Axy+Ayx))/2,
A+0 = (Axz-i Ayz)/sqrt(2),
A00 = Azz,
A-- = (Axx-Ayy+i(Axy+Ayx))/2.
```

All nine formulas are now transcribed in `src/order_parameter_basis.f90` and
covered by known-component and general complex round-trip tests.

There are deliberately two harmonic-to-Cartesian routines. The mathematical
inverse of the nine formulas reconstructs a general complex 3 by 3 field. The
transport reconstruction copied from `new_src/interpol.f90` is a distinct,
source-faithful operation: its formulas enforce

```text
Axy = -Ayx.
```

It is therefore valid only on the restricted axial radial subspace assumed by
that source and must not be called a general inverse. The axial radial
embedding uses this transport reconstruction to reproduce current behavior;
future unconstrained 2D fields remain Cartesian and do not pass through it.

## 5. Fermi-surface and trajectory coordinates

**Source (paper Eqs. 12-13):**

```text
qhat = (q_perp cos(theta_p), q_perp sin(theta_p), q_z),
q_perp = sqrt(1-q_z^2),
n_p = (cos(theta_p),sin(theta_p),0),
s = n_p dot r,
b = e_z dot (r cross n_p).
```

The radial codes sample `theta_p` at azimuthal midpoints with equal weights and
use the 11-point `q_z` quadrature in `gauss11.dat` (weights divided by two on
input). A 2D point and a Fermi-surface direction determine the straight line
used for transport. The transport step along the projected trajectory is not
the same object as a 2D mesh spacing and must become an independent parameter.

## 6. Spin-Nambu ordering and diagonal self-energy

**Source (paper Appendix A):** the Nambu spinor is `Psi=(U,V)^T`, where `U`
and `V` are two-component spinors. For triplet pairing and spin-independent
diagonal self-energy the block matrix is

```text
Sigma_hat = [ nu(qhat),                 (d dot sigma)i sigma_y
              i sigma_y(d* dot sigma), nu(-qhat)               ],
```

with `nu(qhat)=-nu(-qhat)` for the mass-current feedback used here. The source
also allows a spin-vector diagonal self-energy; the current radial calculation
sets that exchange field to zero inside `setenergy1/2`.

The input `F1s` is converted in both legacy and new code to
`F1s/(1+F1s/3)`. Appendix A, Eq. (A10), relates `nu` to `qhat dot j` with an
additional normalization containing the normal-state density. The code arrays
`vx,vy,vz` are the vector coefficients used to construct this diagonal
self-energy/current feedback. They are therefore **current-related mean
fields**, not automatically a physical current density in SI units. A future
observable routine must state the prefactor that converts them to `j`.

## 7. Matsubara and Riccati conventions

**Source:** the current solver reads positive Ozaki poles `zp` and forms
`cen=i*zp`. It integrates two stable Riccati/coherence solutions from opposite
ends of a trajectory. Internally a coherence amplitude is represented by four
complex coefficients `(gamma_0,gamma_x,gamma_y,gamma_z)`, then assembled into
a 2 by 2 spin matrix. Normalization matrices are

```text
N_plus  = (1-gamma gamma_tilde)^(-1),
N_minus = (1-gamma_tilde gamma)^(-1).
```

The reconstructed dimensionless propagator contains explicit factors of `i`.
The source reconstruction is now implemented as explicit 2 by 2 spin blocks
and tested for arbitrary nonsingular coherence matrices and for coherence
values frozen from a source trajectory. Pointwise multiplication gives
`g_hat^2=-1` to roundoff. This fixes the normalization of the compatibility
kernel; it does not by itself select the normalization of future external
interfaces that may use the common `-pi^2` convention.

The interval update in `new_src/riccati.f90:ricc` and its endpoint bulk
coherence calculation are now transcribed as pure, explicit-kind routines.
Their unusual half-step coefficient interpolation is retained exactly. Frozen
forward and reverse checkpoint values are generated by small reference drivers
that link the original module, and the modern routines reproduce every stored
checkpoint. This validates source behavior only; it does not resolve the
general tilde questions above.

The coefficient-to-matrix maps are kept explicitly distinct. The forward
coherence is assembled as
`(gamma_0 I + gamma_alpha sigma_alpha) i sigma_y`. The reverse source branch
uses the exact `Gt` signs from `new_src/riccati.f90:green` and is named as a
legacy reverse conversion rather than being silently declared the general
tilde operation. Tests cover the spin-scalar limit, triplet spin-rotation
covariance, the conjugate-coherence anomalous-block relation, and Nambu
normalization.

**Project proposal:** define the equilibrium tilde operation once at the field
API boundary as

```text
tilde X(qhat,R,z) = X(-qhat,R,-z*)*.
```

For Matsubara `z=i epsilon_n`, this simplifies accordingly. The existing code
uses local complex conjugation when constructing its two propagation branches;
that implementation shortcut must not be promoted to the general tilde
definition without verifying the momentum, energy, transpose, and spin
structure for every Nambu block.

## 8. Self-consistency fields and residuals

The nonlinear state contains 21 real values per spatial point:

```text
18 values = Re/Im of the nine A_{alpha i} components,
 3 values = real vx, vy, vz.
```

The packing order in `packses.f90` is part of the radial regression contract.
The new 2D code should expose typed fields and let an accelerator adapter own
this flattening order.

Both an RMS/average residual and a maximum-component residual are required.
Convergence means that the stated stopping rule is met; reaching `itmax` is a
distinct `iteration_limit` status even when the executable returns zero.

The source `green` contribution has now been extracted behind the explicit
spin-Nambu propagator interface. For one momentum direction and one positive
Ozaki pole it accumulates a complex 3 by 3 gap map and a real three-component
current-related mean field. The source factors `3*awei*pwei*Rp` for the gap and
`0.5*T*awei*pwei*Rp` for the current are preserved. The final source `esum`
factor applies only to the gap map and is therefore a separate finalization
argument. A driver linked to the original module supplies the frozen numerical
fixture.

The historical `green` and `makeprops` declarations marked their accumulating
`tem` argument as `intent(out)` and then read it. That is undefined Fortran
semantics. The working `new_src` declarations have now been corrected to
`intent(inout)` without changing the arithmetic, and the modern accumulator is
likewise initialized explicitly. Its zero-initialized result agrees with the
source routine.

## 9. Density and current observables

**Source (paper Eq. 21):** the pair density used to locate half cores is

```text
n_pair = Tr(A_tilde A_tilde^dagger)/3
       = sum_{alpha,i} |A_tilde_{alpha i}|^2/3.
```

The radial `curr` file writes the square root of the analogous code-basis
quantity. The Python summary plot squares/reconstructs directly from `op_xyz`
and normalizes by its profile maximum. It labels this a **pair-density proxy**.
It must not be called the superfluid density `rho_s`: `rho_s` is a response
coefficient/tensor and requires a separate calculation.

For the 2D double core, physical current should be output as `j_x,j_y`; the
stream function convention in the paper is `j_x=partial_y psi` and
`j_y=-partial_x psi`. Stagnation points satisfy `j=0`.

## 10. Units requiring a freeze decision

The paper uses `T/Tc`, the bulk gap `Delta_0(T)`,
`xi_0=hbar vF/(2 pi Tc)`, and
`R0=(1+F1s/3) xi_0`. The legacy and `new_src` codes clearly use `T/Tc` and a
dimensionless weak-coupling gap, but the source does not state unambiguously in
one place whether printed radii are in `xi_0`, `R0`, or another transport unit.
No new API or plot should attach a coherence-length symbol to the radius until
this is resolved from the original derivation or a known numerical benchmark.

The following must be decided and tested before G1:

- the exact energy unit of `zp`, `Delta`, and `nu`;
- the length unit multiplying the Riccati derivative;
- whether code-basis `A` is paper `A` or `A_tilde` and where `Delta_0` resides;
- normalization conversion at external interfaces (`-1` versus `-pi^2`);
- the full tilde operation for scalar, spin-vector, singlet, and triplet blocks;
- physical prefactors for mass current, free energy, and response-defined
  superfluid density.
