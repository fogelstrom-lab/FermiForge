# Current `new_src` normal-core reference

This directory is the primary radial reference for the first FermiForge 2D
normal-core benchmark. It was generated in an isolated directory with the
current `new_src` executable and was not copied over any working source data.

The calculation used a free unit-winding vortex, `T/Tc=0.30`, `F1s=5.4`, 48
azimuthal directions, Gauss-11 polar quadrature, eight Ozaki poles, and the
100-point tangent grid extending to radius 70. It started from `istart=0`, ran
five simple-relaxation iterations, then converged after 37 Anderson iterations
with `p_max=5` and history limit 10.

The final reported `(average, maximum)` residuals were
`(5.130e-9, 1.562e-6)`, both below the declared `2e-6` tolerance. The complete
solver output and accelerated-iteration history are retained. `manifest.json`
records the executable hash, numerical controls, terminal status, and hashes
of every reference file.

`point_24.expected` is an independent complete `new_src` self-consistency map
at radial index 24. Its comparison with the modern 2D point map is documented
in `../../point_map_comparison.md`.
