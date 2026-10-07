# FermiForge implementation note

## Living hands-on framework guide

`fermiforge_framework_guide.tex` is the standalone LaTeX edition of the
hands-on package map (revision 0.1, 6 October 2026). It includes a vector
flowchart, source-module map, practical exercises and a revision record.
Open it in the native LaTeX editor for an editable source and PDF preview;
no external figures or project includes are required. Update its revision/date
macros and revision table as the implementation evolves. Keep the companion
`../FRAMEWORK_GUIDE.md` navigation guide aligned; synchronization is manual.

## Detailed mathematical note

The canonical PDF is `output/pdf/fermiforge_algorithm_and_architecture.pdf`
relative to the repository root. Revision 0.2 (23 September 2026) includes
the current smooth/extended grid, blended-quadratic ray sampling, evolving
asymptotics, all four seed formulas and a standalone run guide.

Sources:

- `fermiforge_algorithm_and_architecture.tex`: equations and architecture.
- `core_initialization.tex`: exact initializer definitions and distinctions.
- `standalone_run_guide.tex`: software, input inventory, full runnable namelist,
  launch, checkpoint inspection and restart.
- `tools/make_technical_note_figures.py`: reproducible scientific schematics.

From the repository root, with Python plotting dependencies and a TeX Live
installation providing `latexmk`, TikZ and listings:

```sh
python3 tools/make_technical_note_figures.py
cd docs/technical_note
latexmk -pdf -interaction=nonstopmode -halt-on-error fermiforge_algorithm_and_architecture.tex
```

The namelist listing is read directly from `examples/2d_quadrupole_Fs1_0.nml`;
keep the repository layout intact. After inspecting the rendered document,
copy the compiled PDF to the canonical `output/pdf` location. Figures are
illustrations of the implemented formulas, not calculated equilibrium data.

The circular-mask driver is fully 2D despite its historical name. Its required
radial reference remains a compatibility dependency and is explicitly not a
quadrupole accuracy target. No production solver or input changes accompany
this documentation revision.
