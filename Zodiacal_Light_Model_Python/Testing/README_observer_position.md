# ZodiSURF real-observer-position patch

## Why

ZodiSURF's `get_zmod()` places the observer strictly **inside the ecliptic plane** when
computing the zodiacal-light line-of-sight integral:

```python
X0 = Re * COSTHETA
Y0 = Re * SINTHETA
# Z0 is implicitly 0 -- there is no observer out-of-plane term at all
```

`Re` and `Theta` (the observer's Sun-distance and ecliptic longitude) come from a simplified
analytic 2-body Earth orbit, optionally pushed out by a fixed scalar `L2` offset
(`Re * 1.010037`). There is no parameter anywhere in the public API for the observer's
**height above or below the ecliptic** -- i.e. no way to express a real Sun-Earth L2 halo
orbit's Z-axis libration, which is exactly the axis that changes the spacecraft's distance
from the interplanetary dust cloud's own (slightly tilted) symmetry plane.

This patch adds that missing degree of freedom.

## What it does

Adds one new optional parameter, `obs_xyz_au`, to both `get_zmod()` and `zkernel()`: the
observer's true heliocentric position, **ecliptic Cartesian, in AU**, one `[X, Y, Z]` triple
per requested point. When given, it is used directly in place of the internal
`Re`/`Theta`-derived, ecliptic-plane-only position -- Z included. When omitted (the default),
behavior is untouched.

Two files are touched:
- `zodisurf/zkernelpy.py` -- the actual line-of-sight geometry (`zkernel()`)
- `zodisurf/get_zmod.py` -- threads the new parameter through the public `get_zmod()` call

## Install

```bash
git clone https://github.com/rosaliaobrien/skysurf.git
cd skysurf/Zodiacal_Light_Model_Python
pip install -e .
cd zodisurf
patch -p0 < zodisurf_observer_position.diff
```

`patch -p0` must be run from inside the `zodisurf/` directory (the diff's paths are relative
to it, e.g. `--- zkernelpy.py`). You should see:

```
patching file zkernelpy.py
patching file get_zmod.py
```

## Verify

```python
from zodisurf.get_zmod import get_zmod
import numpy as np

# regression check: identical to the unpatched package's own documented test case
print(get_zmod(lambda_=1.25, phase_type="kelsall", day=1.0, lon=0.0, lat=90.0))
# -> [0.11494]

# new capability: an explicit 3D position
print(get_zmod(lambda_=1.076, phase_type="skysurf", day=np.array([12000.0]),
               lon=np.array([0.0]), lat=np.array([90.0]),
               obs_xyz_au=np.array([[1.0, 0.0, 0.003]])))
```

The regression check must reproduce `0.11494` exactly -- the patch is additive only and
changes nothing when `obs_xyz_au` isn't passed.

## Usage

`zodi_analysis.py` calls this automatically when run with `--orbit-file`: it computes each
model-grid time's true heliocentric ecliptic position (real spacecraft ephemeris + Earth's
real position, via `spacecraft_positions()`) and passes it straight through as `obs_xyz_au`
for the ZodiSURF-backed bands (VIS, Y, J). No other script changes are needed once the patch
is applied -- `zodi_analysis.py` degrades gracefully to ZodiSURF's stock `L2=True`
approximation if it isn't.

## Scope of the effect

Checked against a real Euclid ephemeris at the north ecliptic pole: the Z-axis term alone
changes the predicted brightness by **up to ~1.1%, mean ~0.6%** relative to the same position
with Z forced to 0. Small next to a ~30% data/model discrepancy, but real, and worth removing
before quoting that discrepancy precisely.
