# Embedded radial one-point map

This benchmark compares the modern serial 2D point-map path with a complete
`new_src/getnewses.f90:getnewop` calculation at radial index 24
(`r=3.6095394702719363`). Both use 48 azimuths, the 11-point polar table, all
eight Ozaki poles at `T=0.30`, the source energy-dependent Riccati substep
rule, and the same current-feedback coefficient.

The independent legacy values are stored in
`tests/reference/legacy_getnewop_point_24.expected`; their generating driver
and exact build command are in `tests/reference/`.

The modern route first embeds the axial radial state on a square Cartesian
grid with half-width 32, then evaluates the target point. Refining only that
embedding gives:

| Cartesian spacing | max absolute gap error | max relative gap error | max absolute mean-field error |
|---:|---:|---:|---:|
| 0.500 | 6.9064e-4 | 2.2945e-3 | 8.9095e-5 |
| 0.250 | 1.9599e-4 | 6.5113e-4 | 2.5708e-5 |
| 0.125 | 3.6068e-5 | 1.1983e-4 | 6.0659e-6 |

The decreasing error identifies the principal discrepancy as Cartesian
embedding/interpolation error, rather than a spin, quadrature, or component
convention mismatch. The default CTest case uses spacing 0.25 and requires a
relative gap difference below `1e-3` and a mean-field difference below
`5e-5`.

Run a chosen refinement with:

```text
work/full-point-map-build/compare_radial_point_map \
  new_src/op_xyz new_src/curr new_src/gauss11.dat new_src/ozaki.dat \
  512 32
```
