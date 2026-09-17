# FermiForge two-dimensional field sampling

Status: field, axial embedding, straight-ray geometry, source self-energy
projection, legacy Riccati checkpoints, and a complete serial one-point
self-consistency map implemented. Full-grid point batching, MPI collection,
spatial residual diagnostics, an Anderson update, a selectable controlled
free-vortex endpoint, and the first static multiscale rectilinear mesh with a
circular active region are connected. Outer-cutoff/domain convergence and a
converged axial iteration remain.

## Scope

The first 2D slice removes cylindrical symmetry only from spatial field
storage and sampling. It deliberately does not change the Riccati equations,
angular or energy quadrature, nonlinear iteration, or free-vortex asymptotic
condition.

The mesh covers `[x_min,x_max] x [y_min,y_max]` with nodal fields and uses a
deterministic flat point index in which `x` varies fastest. Both uniform
coordinates and explicit nonuniform rectilinear coordinates use the same
field and interpolation interfaces. A circular active-region mask may limit
self-consistency mapping and residuals while leaving a frozen exterior halo
available to every trajectory. The typed state stores

```text
A(spin, orbital, point)       complex, 3 x 3 x npoint
nu(component, point)          real,    3 x npoint
```

`nu` denotes the current-related mean field represented by `vx,vy,vz` in the
radial source. It is not labelled as physical current or superfluid density.

## Interpolation contract

The initial backend is local bilinear (`Q1`) interpolation. A stencil contains
four flat node indices and four real weights. A point on any domain edge is
inside; a point outside is reported explicitly rather than extrapolated. Only
round-off-sized excursions are clamped to the edge.

Complex Cartesian components are interpolated directly. Amplitude and phase
are not separated. Cell lookup is by binary search in each explicit coordinate
array, and interpolation uses the local cell widths. Constant, affine, and
bilinear manufactured fields must be reproduced to roundoff; smooth quadratic
and axial-vortex fields must converge under refinement.

The stencil is a public value type so fixed-mesh ray geometry can precompute
indices and weights. Nonlinear iterations then load fields and apply weights
without repeating cell searches.

## Iteration-vector contract

The typed state does not expose a solver-dependent memory layout. The first
adapter explicitly reproduces `new_src/packses.f90`: component-major blocks of
the real and imaginary parts of the nine order-parameter components followed
by three real mean-field blocks. Its length is `21*npoint`. Unlike the older
fixed-form iterator contract, this `new_src` layout does not multiply the mean
fields by ten.

## Translation-invariant trajectory contract

For a target point `r0=(x0,y0)` and Fermi-surface momentum
`p=(px,py,pz)`, a straight projected trajectory is

```text
r_xy(s) = r0 + s*(px,py).
```

The spatial `z` coordinate is absent. The momentum component `pz` is retained
and the sampled triplet pair potential is

```text
Delta_alpha(p,r) = sum_i A(alpha,i,r) p_i.
```

The straight-ray layer now intersects this line with a rectangular 2D domain,
places the target exactly at `s=0`, limits every interval to a caller-provided
maximum step, and caches one bilinear stencil per path sample. A direction
parallel to the invariant `z` axis has one spatial sample. The full normalized
three-component momentum is stored even though only `(px,py)` controls the
projected geometry.

Straight-ray geometry remains independent of boundary physics. The first
free-vortex endpoint policy extends its sampled self-energy path beyond the
rectangle without adding fictitious interpolation stencils; later wall
reflection remains a separate boundary and transport policy.

## Free-vortex outer continuation

The selectable free-vortex endpoint uses a rectangular numerical field domain
and an enclosing outer circle. At every exterior point it projects radially
back to the rectangle, samples the numerical boundary field, and reconstructs
an analytic tail. The bulk part is

```text
A_bulk = Delta_bulk exp(i N phi) I.
```

Departures with exactly one `z` index decay as `r_match/r`; all other
order-parameter departures decay as `(r_match/r)^2`. The current-related mean
field decays as `r_match/r`. These are the explicit Cartesian powers used by
the axial free-vortex extrapolation in `new_src/interpol.f90`. The construction
is continuous at the rectangular boundary and allows a fully 2D boundary
field; it does not impose cylindrical symmetry inside the box.

The path is extended with a separately bounded tail step. The two stable
Riccati branches are initialized on the outer circle and propagated through
the analytic tail and numerical field to the target. The previous local
endpoint remains selectable for regression. A uniform bulk-B field is
unchanged by the continuation, the tail powers and boundary continuity have
unit tests, and the MPI result agrees with the serial result for 1, 4, and 10
ranks. Production acceptance still requires refinement in box size, outer
radius, and both interior and tail trajectory steps.

## Trajectory self-energy contract

For each cached path sample the current compatibility layer constructs

```text
d_alpha(s) = sum_i A(alpha,i,r(s)) p_i,
nu(s)      = a_1 dot(v_current(r(s)),p).
```

Here `a_1` is the caller-supplied current-feedback scale corresponding to the
transformed `aa0` used by `new_src`; no physical-current interpretation is
attached to the stored vector. The two legacy propagation branches place `d`
and its complex conjugate in the same coefficient slots as `setenergy1` and
`setenergy2`, use `z-nu`, and explicitly retain a zero spin-exchange field.

The side-effect-free interval kernel is an arithmetic-preserving transcription
of `new_src/riccati.f90:ricc`. Variable 2D path intervals are allowed, while the
RK4 substep count is a separate control. A trajectory driver initializes each
end with the source bulk coherence and propagates both stable branches to the
target. Optional boundary-relaxation distance is explicit rather than hidden
as `10*dx`.

## Radial reference embedding

An adapter converts an axial `new_src` radial profile to harmonics, applies the
source phase `exp(i*(M-sigma-ell)*phi)`, reproduces the source three-point
radial interpolation including negative-radius parity at the origin, and
embeds the result on Cartesian nodes. It reports a Cartesian node beyond the
radial support instead of extrapolating silently.

The adapter intentionally uses the restricted reconstruction copied from
`new_src/interpol.f90`, which enforces `Axy=-Ayx`. A separate exact inverse
exists for general basis conversions. Fully 2D fields are stored and sampled
directly in Cartesian form, so the axial restriction is not imposed on the
new solver state.

## Validation and integration order

1. Run the manufactured mesh and interpolation tests.
2. Freeze the radial harmonic-to-Cartesian convention and embed a trusted
   `new_src` radial solution on the 2D nodes. (Implemented and unit tested;
   numerical-oracle comparison remains part of step 3.)
3. Compare sampled self-energies along selected cached rays with radial
   interpolation. (Manufactured projection and radial-embedding refinement
   tests pass.)
4. Extract a side-effect-free Riccati kernel and compare propagation
   checkpoints with the source routine. (Complete for the legacy coefficient
   kernel in both propagation directions.)
5. Compare one complete 2D self-consistency-map evaluation with the radial
   map before enabling 2D relaxation or a double-core seed. (Complete at
   radial index 24: the 528-direction, eight-pole map has `6.51e-4` relative
   error at Cartesian spacing 0.25 and decreases to `1.20e-4` at spacing
   0.125.)

The implemented multiscale bridge is a fixed, nonuniform tensor-product mesh
with user-selected fine, medium, and outer spacings. It validates nonuniform
sampling and separates an active circular solve region from a frozen halo,
but it is not error-driven AMR. See `MULTISCALE_RECTILINEAR_MESH.md`.

Adaptive block meshes remain an outer solve-estimate-refine-transfer
operation. A mesh change invalidates cached ray stencils and nonlinear-history
vectors; both must be rebuilt or reset.
