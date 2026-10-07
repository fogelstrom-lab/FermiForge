# Importing legacy two-dimensional fields

The old two-dimensional solver writes one row of the complex order-parameter
matrix to each of `op_x`, `op_y`, and `op_z`. Every record contains `x`, `y`,
and three real/imaginary pairs. The `curr` records contain

```text
x y gap_norm in_plane_magnitude mean_field_x mean_field_y mean_field_z
```

FermiForge reads these four primary files directly. Blank row separators are
optional, coordinates must describe an x-fast rectilinear grid, and all four
files must contain identical coordinates. The importer checks the stored gap
norm and in-plane magnitude while preserving the three mean-field components.

The raw `2D_benchmarks/` directory is deliberately ignored by Git because it
is archived separately and is hundreds of megabytes. With a local copy in the
project directory, configure and build FermiForge, then inspect it with:

```text
work/modern-build/import_legacy_split_field_map_2d \
  '2D_benchmarks/Double_core_vortex_T=0.30_Fs1=5.4'
```

To convert the state into the versioned FermiForge field-map format, supply an
output path under the ignored `work/` directory:

```text
work/modern-build/import_legacy_split_field_map_2d \
  '2D_benchmarks/Double_core_vortex_T=0.30_Fs1=5.4' \
  work/double_core_fields_2d.dat
```

## Generate the comparison figures

The archived double-core state can be plotted directly, without first writing
the large intermediate field map:

```text
tools/plot_legacy_double_core_reference.py \
  '2D_benchmarks/Double_core_vortex_T=0.30_Fs1=5.4'
```

The command writes twelve PNG files and a machine-readable manifest under
`2D_benchmarks/Double_core_vortex_T=0.30_Fs1=5.4/plots/`. Five full-cell
figures correspond to the modern final-state Cartesian amplitude, Cartesian
phase cosine, orthonormal harmonic amplitude, harmonic phase cosine, and pair
density/current-related-mean-field figures. The same five figures are also
written for the default `[-30,30]^2 xi0` core window, so that the double core
is not visually compressed by the full `[-60,60]^2 xi0` legacy domain. The
remaining figures show the dcvlong-Fig.-4 symmetry-axis profiles and the
preliminary/Anderson convergence records from the final block of `qcv.log`.

The plotting tool verifies all coordinates, the redundant gap norm and
in-plane field magnitude, and the normalized `dens` quantity. For this archive
the two hard-core minima are at `y=+/-11.8 xi0`, giving a separation of
`23.6 xi0`. The vector in `curr` is labelled as a current-related Fermi-liquid
mean field; it is not silently identified with a physically normalized mass
current. There is no archived initial field or iteration-resolved field
history, so initial-state, structure-history, and modern probe/shadow figures
cannot be reconstructed from this final state.

The auxiliary `Cpp` through `Cmm` files are diagnostics rather than restart
state. They use the historical two-dimensional postprocessor convention. The
five unmixed correspondences are direct:

```text
FermiForge:  A++  A+-  A00  A-+  A--
legacy:      Cpp  Cpm  Coo  Cmp  Cmm
```

For the four mixed-zero channels, both the filename ordering and amplitude
normalization differ:

```text
FermiForge amplitude:  |A+0|  |A0+|  |A0-|  |A-0|
legacy amplitude:      Cmo/2  Com/2  Cop/2  Cpo/2
```

The corresponding legacy cosine-of-phase column uses the same permuted file
without the factor of two. These files must therefore not be read as
same-named FermiForge harmonic components. The reference plotting tool avoids
the ambiguity by recomputing all harmonics from `op_x`, `op_y`, and `op_z`.
