# Continuing with full tensors and evolving asymptotics

Run from the project root. This preserves the saved interior field and resets
Anderson history. The halo is regenerated from the current state before the
first map and after each update. Both incoming trajectory tails and halo use
the common free-vortex asymptotic implementation, without a core-type switch.

```bash
caffeinate -i env OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,sm \
  python3 tools/run_axisymmetric_core_benchmark.py \
  --input runs/260922-152321-a-phase-continued-radial-reference/input.nml \
  --restart runs/260922-152321-a-phase-continued-radial-reference/final_fields_2d.dat \
  --endpoint-policy state_asymptotic \
  --fit-inner-radius 10 --matching-radius 14 \
  --build-dir work/double-core-scratch-build \
  --case-name a-phase-evolving-asymptotics \
  --ranks 10 --max-iterations 5 --checkpoint-interval 1 --formats png
```

The matching circle at 14 lies inside the active radius 16, leaving room for
interpolation. These radii are numerical test settings, not core-type constants;
their adequacy still needs a matching-radius/domain convergence study.
The bulk reference remains the specified phase-wound bulk B state. Departures
in the mixed z/in-plane block decay as 1/r; other blocks as 1/r^2. The Fermi-liquid
field uses both powers fitted at radii 10 and 14. Coefficients depend on angle
and the current field. This is an asymptotic closure, not an independent solve
of the exterior region.

The physical radial embedding now preserves all nine components and rotates
them as a tensor. Archived fields are not silently rewritten: this restart
still contains the old evolved interior, which may initially readjust.
New comparison plots use metrics.txt.radial_reference.dat, freshly embedded
with the corrected transform. Old radial transport comparison tests explicitly
retain the historical Axy=-Ayx projection; they test legacy reproduction only.

Plots are attempted from saved results even after a launcher error, while the
runner retains a nonzero exit status and records launcher_error in its manifest.
This does not fix the outstanding MPI shutdown warning. Do not infer physical
convergence from successful plotting or from benchmark passed alone.
