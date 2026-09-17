# Legacy checkpoint generators

These small programs call the original `new_src/riccati.f90` routines and are
the provenance for frozen values in the modern Riccati tests. They are not part
of the normal CMake test build because they intentionally depend on the legacy
module, global state, and default-kind promotion flags.

Regenerate the interval checkpoint with:

```text
/opt/homebrew/bin/gfortran -std=f2018 \
  -fdefault-real-8 -fdefault-double-8 -fimplicit-none \
  -Inew_src tests/reference/legacy_riccati_checkpoint_driver.f90 \
  new_src/riccati.o new_src/global_dec.o \
  -o /tmp/fermiforge_legacy_riccati_checkpoint
/tmp/fermiforge_legacy_riccati_checkpoint
```

Replace the driver filename/output name with
`legacy_riccati_trajectory_driver.f90` and
`fermiforge_legacy_riccati_trajectory` for the five-point forward/reverse
trajectory. Frozen values must be updated only when a source behavior change
has been reviewed and accepted, never merely to make a regression pass.

`legacy_green_integrand_driver.f90` links the original `new_src` modules and
freezes one contribution from `riccati.f90:green`. The legacy routine declares
its accumulating array as `intent(inout)`. This is the corrected declaration
for an array that the source always used as an accumulator; the earlier
`intent(out)` declaration was invalid Fortran. The modern integrand is compared
with the zero-initialized legacy result.

Regenerate that fixture with:

```text
mkdir -p work/legacy-green-reference
PATH=/opt/homebrew/bin:/usr/bin:/bin /opt/homebrew/bin/mpifort \
  -std=f2018 -fdefault-real-8 -fdefault-double-8 -fimplicit-none \
  -Jwork/legacy-green-reference -Iwork/legacy-green-reference \
  new_src/global_dec.f90 new_src/riccati.f90 \
  tests/reference/legacy_green_integrand_driver.f90 \
  -o work/legacy-green-reference/legacy_green_integrand
work/legacy-green-reference/legacy_green_integrand
```

`legacy_getnewop_point_driver.f90` runs the complete original angular and
Ozaki loop at one radial grid point, without MPI reduction or nonlinear
mixing. It is the independent reference for the upcoming 2D point-map
comparison. It defaults to radial index 24; pass another index as its first
argument if required. Build it with the source's required default-kind
promotion, then run it in a scratch directory containing copies of the input
tables and radial state:

```text
mkdir -p work/legacy-getnewop-reference/build \
         work/legacy-getnewop-reference/run
PATH=/opt/homebrew/bin:/usr/bin:/bin /opt/homebrew/bin/mpifort \
  -std=f2018 -fdefault-real-8 -fdefault-double-8 -fimplicit-none -O0 \
  -Jwork/legacy-getnewop-reference/build \
  -Iwork/legacy-getnewop-reference/build \
  new_src/global_dec.f90 new_src/bulkgap.f90 new_src/init_calc.f90 \
  new_src/interpol.f90 new_src/riccati.f90 \
  tests/reference/legacy_getnewop_point_driver.f90 \
  -o work/legacy-getnewop-reference/build/legacy_getnewop_point
cp new_src/qcv.inp new_src/gauss11.dat new_src/ozaki.dat \
   new_src/op_xyz new_src/curr work/legacy-getnewop-reference/run/
cd work/legacy-getnewop-reference/run
../build/legacy_getnewop_point 24 < qcv.inp | tee reference-point-24.txt
```
