# Free-form cylindrical trajectory port

`trajectories.f90` defines the `cylindrical_trajectories` module, with public
`trajectories` and private `zigzag`. `poretraj.f90` ports the interactive test
program from `poretraj.f` and uses that module instead of embedding another
copy of the transport geometry. The original `.f` sources are unchanged.

This is a behavior-preserving geometry port. The active `new_src` radial
solver now uses this module for `icyl=1`; see
[`RUNNING_RADIAL_SPECULAR_CYLINDER.md`](../../docs/RUNNING_RADIAL_SPECULAR_CYLINDER.md). The modern full 2D solver does
not yet use this boundary path.
It models specular reflections at the side wall of an infinite circular
cylinder; there are no end caps or order-parameter boundary conditions here.

## Build and run

From the repository root, build in an isolated directory:

```sh
mkdir -p work/cylinder-demo
cd work/cylinder-demo
gfortran -std=f2018 -fimplicit-none -fcheck=all \
  ../../new_src/3DFS_MPICodes/trajectories.f90 \
  ../../new_src/3DFS_MPICodes/poretraj.f90 -o poretraj
./poretraj
```

The program prompts for radius, integration step, number of steps, initial
position, and unit momentum direction. It writes the original `traj` format
in the current directory. Use a strictly interior initial position.

## Compatibility and limitations

- The argument ordering matches `trajectories.f`, with arrays indexed
  `0:nsn` instead of the fixed `0:20000` limit. The driver allocates its arrays
  dynamically and eliminates the original COMMON blocks.
- Default `real` is deliberate: compile with `-fdefault-real-8
  -fdefault-double-8` when matching the existing radial solver. A later
  explicit-kind API should have its own regression gate.
- The original reflection arithmetic, branch choices, and sampled momentum
  sign conventions are retained. In particular, exactly transverse rays
  change the caller's `pz` to `+/-1e-10`; coordinate increments use the original
  direction. The public interface explicitly marks `pz` as `intent(inout)`.
- The legacy near-axis epsilon substitution and wall/grazing degeneracies
  remain. Exact wall starts are excluded by the driver; the library expects
  valid inputs. Long reflection sequences accumulate rounding error.
- The legacy Makefile and the production CMake solver targets are unchanged.
  The active `new_src/Makefile` and isolated runner now build this module.
  Compile the module before the new driver; do not link both legacy and
  modern drivers into the same executable.

## Regression

Run the short standalone check with GNU Fortran available:

```sh
python3 tests/check_cylindrical_trajectories.py
```

The check builds in a temporary directory with bounds checks and floating-point
exception traps, in both default-real and real-8 modes. Nine cases compare all
12 sampled coordinate/momentum arrays and the `pz` side effect exactly against
`trajectories.f`: oblique rays in both axial directions, transverse rays,
axial rays, exact diameter wall hits, multiple reflections per step, and a
near-wall start. Finite values, cylinder containment, and momentum norms are
also checked with precision-dependent tolerances. Four driver cases in each
precision mode compare the complete `traj` output byte for byte with the
original standalone `poretraj.f`.
