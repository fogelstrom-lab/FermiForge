# Running and plotting the FermiForge 2D demonstration

## What this test is

The current demonstration tests the Cartesian mesh, bilinear interpolation,
field storage, output format, and plotting workflow. It generates a smooth
manufactured singly quantized vortex with a circulating manufactured current
and a weak axial current.

It is **not** a self-consistent quasiclassical solution. The trusted
`new_src` self-consistency calculation has not yet been connected to the 2D
mesh.

## Base directory

Run commands from

```text
/Users/mikael/Documents/Codex/3he-vortex-modernization
```

## Simplest complete test

```text
/opt/homebrew/bin/python3 tools/run_2d_demo.py
```

The script uses the normally selected Apple developer toolchain, builds the
Fortran driver, creates a unique timestamped directory under `runs/`, archives
the input and logs, writes the 2D field map, and produces PNG and PDF figures.

The terminal prints the exact run directory. It contains:

```text
input.nml
fields_2d.dat
manifest.json
run.log
plot.log
plots/
```

## Making a test input

Copy and edit the supplied namelist:

```text
cp examples/2d_sampler_demo.nml work/my_2d_test.nml
```

Then run it with a descriptive case name:

```text
/opt/homebrew/bin/python3 tools/run_2d_demo.py \
  --input work/my_2d_test.nml \
  --case-name fine-grid
```

The input parameters are:

| Name | Meaning |
|---|---|
| `x_minimum`, `x_maximum` | horizontal domain bounds |
| `number_of_x_cells` | number of uniform horizontal cells |
| `y_minimum`, `y_maximum` | vertical domain bounds |
| `number_of_y_cells` | number of uniform vertical cells |
| `core_width` | manufactured vortex-core scale |
| `texture_width` | radius of the manufactured off-diagonal texture |
| `off_diagonal_scale` | strength of off-diagonal order-parameter components |
| `current_scale` | circulating in-plane manufactured current |
| `axial_current_scale` | weak manufactured current along the invariant axis |
| `phase_winding` | integer phase winding |
| `output_at_cell_centers` | if true, exercise interpolation at cell centers |

The `output_file` entry is used when the Fortran program is run directly. The
Python wrapper overrides it with a collision-free path in the run archive.

## Plotting an existing map

```text
/opt/homebrew/bin/python3 tools/plot_2d_fields.py \
  runs/RUN_NAME/fields_2d.dat
```

Use `--output-dir`, `--prefix`, `--formats`, `--phase-mask-fraction`, and
`--quiver-target` to adjust output names and display density.

Three figures are produced:

1. A 3 by 3 grid of `|A_mn|`. All panels use one shared scale. Zero is white;
   increasing amplitude progresses through a monotonic light-purple scale to
   dark indigo.
2. A 3 by 3 grid of `cos(arg(A_mn))`. Negative values are blue, zero is white,
   and positive values are red. Grey means the amplitude is too small for its
   phase to be meaningful.
3. Pair density, in-plane current magnitude with direction arrows, and a
   separate signed axial-current panel. The axial panel has its own symmetric
   scale so a weak `j_z` remains visible.

The blue-white-red scale avoids a red-green distinction, the amplitude scale
has monotonic lightness, and every panel includes numerical colorbar ticks.

## Running only the compiled Fortran program

```text
cmake -S . -B work/fermiforge-2d-build -G Ninja \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran
cmake --build work/fermiforge-2d-build --parallel 10 \
  --target fermiforge_2d_demo
./work/fermiforge-2d-build/fermiforge_2d_demo \
  examples/2d_sampler_demo.nml work/fields_2d.dat
```

The ASCII output declares its columns in the header. The same format is the
target for the later self-consistent 2D solver, allowing the plotting tool to
remain unchanged when manufactured fields are replaced by physical results.
