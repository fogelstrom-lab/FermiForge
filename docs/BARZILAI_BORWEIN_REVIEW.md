# Review of new_src/iter_BB.f90 (2026-09-27)

Update: a separate modern implementation is now available; see
[BB iteration](BB_ITERATION.md). The review below describes the original source.

This is a source-code review, not a validated 2D BB implementation.
The original routine and its alternative iter_BB_ok.f90 have not been changed.

## What the routine actually computes

The real vector includes all real/imaginary order-parameter components and
the three real Fermi-liquid fields. packses.f90, mode 2, explicitly produces
the residual r = F(X) - X. iterateBB updates X_next = X + alpha*r.
The first update is an undamped substitution, alpha=1.

Writing d = r_previous - r_current, later updates alternate the multipliers

- ||r_previous||² / |d dot r_previous|,
- |d dot r_previous| / ||d||².

The multiplier is applied to the previous alpha, then alpha is clipped to
[0.001,100]. The second iteration uses the second expression; subsequent
iterations alternate. When the preceding displacement really equals
alpha_previous*r_previous, these expressions are the BB1/BB2 secant step
lengths for the root residual X-F(X), with absolute curvature substituted
for signed curvature. This equivalence follows by substituting that
displacement into the two secant quotients. It need not hold after projection,
clipping of individual components, remeshing, or switching iteration engines.

## Issues before reuse

1. alpha is declared intent(out), but is read in alpha=alpha*delta on later
   calls without first being defined in that call. Caller storage may happen
   to retain its bits, but Fortran does not guarantee this value. A modern
   implementation should own the previous step in explicit iterator state.
2. Zero or tiny ||d||² and |d dot r_previous| are unguarded. Clipping after
   division is not a reliable NaN/overflow safeguard.
3. The convergence normalization divides by ||X|| without a zero guard.
   tolA is a relative Euclidean norm; tolL is N*max|r|/||X||, not a raw
   average and maximum error. For comparison use the modern common diagnostics.
4. Absolute curvature forces positive steps even for negative secant curvature.
   Preserve that choice in a clearly labeled historical mode if benchmarking
   it; do not silently present it as identical to a safeguarded signed method.
5. alpha=100 can produce very large moves. There is no residual-growth
   acceptance/rejection check. Its initial full step is also much more
   aggressive than the AA startup factor 0.01.
6. Saved allocatable buffers are allocated whenever iiter=1, without a
   deallocation/reset path. This is unsuitable for repeated in-process restarts.

## Assessment and proposed next experiment

BB is a credible separate option: it needs a small amount of vector history,
dot products and no extra map evaluation merely to propose a step. The
scalar step can grow when the residual varies slowly along the update
direction. That makes it relevant to the observed slow core expansion, but
does not establish that it will accelerate this particular mode or converge.

First isolate AA history resets using the same starting field. Then implement
BB as an independent modern engine, with actual successive active-vector
displacements, explicit state/reset, finite/tiny-denominator checks,
configurable step bounds, and logged fallback decisions. Include both OP
and FL components in the same packing and retain existing halo refreshes.
Reset BB history after a grid change or an engine switch.

Compare from the same checkpoint, with unchanged grid, quadrature, asymptotic
policy and residual definition. Record core-size evolution as well as RMS,
maximum and spatial residuals. Residual reduction alone does not show that
the slow physical core-size mode has equilibrated.

## Observed AA20/simple3 run

In runs/260926-205451-single-harmonic-0plus-cycled-aa20-simple3:
RMS residual is 1.4360e-4 at update 24 and 1.3349e-5 at update 30,
a factor of about 10.8. Across the entire run it goes from 2.7682e-5 to
9.2392e-6 (about 3-fold); it is not monotonic. The residual at update 21
is measured before the first simple step, so its spike cannot be attributed
to that as-yet unapplied simple update.
