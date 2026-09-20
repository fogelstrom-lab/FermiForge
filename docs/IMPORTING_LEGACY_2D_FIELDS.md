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

The auxiliary `Cpp` through `Cmm` files are diagnostics rather than restart
state. They use the historical two-dimensional postprocessor convention: the
mixed zero/plus/minus amplitudes differ by a factor of two from the
orthonormal spherical basis used by FermiForge and `new_src/op_harm`. They
must therefore be compared through an explicitly labelled legacy convention;
they must not be read as FermiForge harmonic components.
