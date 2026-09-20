# Historical `nop`, `aop`, and `dop` vortex seeds

The authoritative source is
`new_src/how_I_initialise_the_different_cores.f`. FermiForge retains that file
unchanged and transcribes its three formulas in
`src/historical_core_seed_2d.f90`. The modern implementation does not infer
missing components: it explicitly zeros the full 3 by 3 order parameter and
all three Fermi-liquid mean fields before applying each formula.

Let `Delta=deltat`, `t=T/Tc`, `aa0=F1s/(1+F1s/3)`, `phi=atan2(y,x)`, and

```text
f(s) = tanh(s)/s, with f(0)=1.
```

The historical seed names are selected through `initialization_mode`.

## `historical_nop`

This is the normal-core seed. With

```text
xgl = 1/sqrt(1-t^2),  s = r/(3*xgl),
```

the only nonzero components are

```text
A_xx = A_yy = A_zz = Delta exp(i phi) tanh(s).
```

The order parameter vanishes at the origin.

## `historical_aop`

This uses the same diagonal vortex and radial scale as `nop`, plus

```text
A_zx =  Delta f(s) [(cos(2phi)+1) + i sin(2phi)]/2,
A_xz = -A_zx,
A_zy =  Delta f(s) [sin(2phi) + i(1-cos(2phi))]/2,
A_yz = -A_zy.
```

The formula seeds the rotationally symmetric A-phase-core branch.

## `historical_dop`

Only this seed broadens its radial scale with the feedback parameter:

```text
xgl = (1+aa0)/sqrt(1-t^2),  s = r/(3*xgl).
```

It uses the same diagonal vortex but adds only

```text
A_zx =  Delta f(s) cos(phi)^2,
A_xz = -A_zx.
```

This is a centered seed with an explicit twofold perturbation. It does not
place two half cores at prescribed coordinates. Their formation and separation
are outcomes of the nonlinear iteration.

The origin convention is also retained exactly: `phi=0`, the diagonal vortex
vanishes, `f(0)=1`, and therefore `A_zx=-A_xz=Delta` for `aop` and `dop`.

## Relation to the London double-core seed

`regularized_london` is a separate modern initializer. It explicitly places
two regularized phase singularities and a finite-width wall between them. It
is useful for continuing an already split branch, but it is not a replacement
for the historical `dop` symmetry-breaking experiment.

The large historical `dop` calculation is launched with:

```text
caffeinate -i tools/run_historical_dop_overnight.py
```

It uses the same cell, mesh, active disk, quadrature, endpoint continuation,
Anderson settings, and checkpoint policy as the London-seeded overnight run.
Thus the initial seed is the controlled difference between the two runs.

For `nop` or `aop`, copy the supplied historical-dop input and change only:

```text
initialization_mode = 'historical_nop'
```

or

```text
initialization_mode = 'historical_aop'
```

`seed_reduced_temperature` must agree with the temperature represented by the
selected Ozaki table. The present supplied tables and examples use `0.30`.
