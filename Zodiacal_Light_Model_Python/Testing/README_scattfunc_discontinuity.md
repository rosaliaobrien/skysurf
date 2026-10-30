# Small, spurious brightness jumps for an out-of-ecliptic observer

**Patch:** `zodisurf_scattfunc_discontinuity.diff`
**Affects:** `zodisurf/zkernelpy.py` — `simpint()` and the `itest` branch in `scattfunc()`
**Severity:** small (order 1e-4 to 1e-3 relative, i.e. roughly 0.01-0.1% of the flux, in our tests) but real and reproducible. Does not affect the stock, in-ecliptic-plane (`L2=True`, or the original 2D approximation) usage at all.

## Symptom

Evaluating ZodiSURF on a fine, regular time grid for a **fixed sky position but a moving,
out-of-ecliptic observer** (e.g. a real spacecraft ephemeris fed into `zkernel()`'s
`obs_xyz_au` parameter, rather than the built-in flat L2 approximation) produces a smooth
curve almost everywhere, punctuated by isolated single-time-step jumps of order
1e-4 to 1e-3 relative. In our case (evaluating a fixed VIS bandpass toward a
high-ecliptic-latitude field, once a week over 3 years, for a spacecraft in a halo orbit
near Sun-Earth L2) these jumps recurred roughly every 5-7 months.

We first suspected a plotting artifact (too coarse a time grid, straight lines between
points), then a spline artifact in our own interpolation of the observer's position. Both
were ruled out (see "How we diagnosed this" below); the effect is intrinsic to ZodiSURF's
own computation and appears even at hourly time resolution.

## Root cause

1. **`scattfunc()`** picks between two geometrically-equivalent phase-angle solutions along
   each line-of-sight integration point with a hard switch:
   ```python
   itest = (LOS >= Re * np.cos(SolElong))
   scat_ang = np.where(itest, np.pi - phase_ang, phase_ang)
   ```
   `scat_ang` itself is continuous across this boundary (`phase_ang == pi - phase_ang`
   exactly at `LOS == Re*cos(SolElong)`), but its *derivative* with respect to `LOS` is
   not — the scattering-angle profile along the line of sight has a genuine slope kink
   at this point.

2. **`simpint()`** integrates along the line of sight with a *fixed step size*
   (`stepsize=0.025` AU by default) out to a bound `b` (called `RANGE` in `zkernel()`).
   The number and absolute placement of the quadrature nodes therefore depend on `b`
   alone. `b` derives from `Re_true = sqrt(X0**2 + Y0**2 + Z0**2)`, the observer's own
   heliocentric distance — so as the observer moves, `b` drifts smoothly, and the discrete
   quadrature nodes drift with it.

3. Whenever a drifting node crosses the (independently-computed, effectively fixed)
   kink location from step 1, the Simpson-weighted sum picks up or drops that node's
   contribution from the "wrong" side of the kink, producing a small, real jump in the
   integrated brightness for that one time step, then it's smooth again.

This is invisible for the stock flat/in-ecliptic-plane observer, because there `Re_true`
only tracks Earth's own slow (~1 year) orbital-eccentricity variation, and in our tests
never happened to sweep a node across the kink at all. Feeding in a real, out-of-ecliptic
observer position (e.g. a halo orbit, which librates in Z on a much shorter timescale,
~180 days for a Sun-Earth L2 halo orbit) makes `Re_true`, and hence the quadrature grid,
drift much faster and in both directions — so the crossing happens often enough to be
clearly visible in a multi-year time series.

## How we diagnosed this

1. Evaluated the model at 6-hour and 1-hour time resolution around several of the jumps.
   The jump persisted at the same timestamp regardless of resolution, and its shape did
   not resemble a coarse-grid plotting artifact.
2. Repeated the same evaluation with `obs_xyz_au=None` (flat L2 approximation). The jumps
   disappeared essentially completely (curvature dropped by 2-3 orders of magnitude).
3. Instrumented `zkernel()` directly (temporary print statements, not included in this
   patch) to dump every density component — `zcloud`, all four `migband` calls, and
   `solring` — at the two time steps straddling a jump. All of them were smooth to within
   their normal step-to-step variation. Only `Scatt` (`scattfunc()`'s return value) showed
   a clean step.
4. Reduced `simpint`'s `stepsize` from 0.025 to 0.005 AU (5x finer). The jump shrank by a
   comparable factor (~4.7x) at the *same* timestamp — the standard signature of a
   quadrature-discretization artifact, not a physical feature (which would not shrink when
   integrated more finely).

## What this patch does — and does not — do

This patch is **documentation only**. It adds comments at the two relevant locations
(reproduced above) explaining the issue, for anyone who next touches this code or hits the
same symptom. **It does not change any numerical output.** We deliberately did not change
`simpint`'s default `stepsize`, since its docstring states it is chosen to
"match the IDL implementation exactly including precision" — we did not want to silently
break that intentional parity for existing users without your input.

## Possible fixes, for your consideration

We see a few options, in roughly increasing order of invasiveness:

- **Finer default `stepsize`.** Simplest, but changes the reference/legacy numerical
  behaviour for *all* callers, including the original flat-observer use case, and costs
  more compute. Our test (0.025 -> 0.005 AU) reduced but did not eliminate the jump.
- **Explicitly insert a quadrature node at the kink.** `LOS = Re*cos(SolElong)` is known
  in closed form before `simpint` is called, so a node could be placed exactly there (or
  the interval could be split into two sub-integrals on either side of it). This would
  fix the issue at its source without changing the node density elsewhere, but touches
  the integration scheme more directly.
- **Replace the hard `np.where` branch with a smooth blend** of the two `scat_ang`
  solutions over a small window around the switch. Removes the kink itself rather than
  compensating for it, but changes results at the sub-quadrature-node level everywhere,
  not just near the switch.
- **Leave as is, documented.** The effect is small (single time steps, ~0.01-0.1% of the
  flux) and only matters for callers feeding in real out-of-ecliptic observer positions,
  which is not the original use case. A code comment may be all that's warranted.

We are not in a position to judge which of these best fits your validation requirements
against the original IDL code, so we've left the choice — and any implementation — to you.

## Reproducing it

```python
import numpy as np
from astropy.time import Time
from zodisurf.get_zmod import get_zmod

# any fixed sky position with a nonzero ecliptic latitude; any real, out-of-ecliptic-plane
# obs_xyz_au(t) with an out-of-plane (Z) excursion that changes sign one or more times
# over the requested time span. We used a Euclid-like halo-orbit ephemeris near Sun-Earth
# L2; a hand-built toy trajectory with a sinusoidal Z(t) should reproduce it just as well.

grid = Time(np.arange(t0, t1, 1 / 24), format="mjd", scale="utc")   # hourly, spanning >= 1 halo period
obs_xyz_au = ...   # shape (N, 3), heliocentric ecliptic AU, with Z crossing zero at least once
val = get_zmod(lambda_=np.full(grid.size, 0.7), phase_type="skysurf",
                day=(grid.jd - Time("1990-01-01").jd + 1.0), lon=..., lat=..., L2=False,
                obs_xyz_au=obs_xyz_au)

d2 = np.diff(val, n=2)
print(np.max(np.abs(d2)) / np.median(np.abs(d2)))   # >> 1 at one or more isolated points
```

## Contact

Happy to share our full test scripts or discuss further if useful.
