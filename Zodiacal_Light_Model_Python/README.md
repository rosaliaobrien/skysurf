# ZodiSURF: SKYSURF Zodiacal Light Modeling Package

Python implementation of the Kelsall et al. (1998) zodiacal light model with SKYSURF enhancements (O'Brien et al. 2025).

## Features

- Ability to run the Kelsall et al. (1998) version and the SKYSURF version (O'Brien+2025).
- **Wavelength Range**: 1.25-240 μm for the Kelsall model; 0.2-1.7 μm for the SKYSURF model.
- **High Accuracy**: 0.0002% difference from IDL package.
- **Vectorized Processing**: Scalar and array inputs

## Installation

```bash
git clone <repository-url>
cd ZodiModel
pip install -e .
```

**Requirements:** Python ≥3.7, NumPy ≥1.19, SciPy ≥1.7

## Usage

```python
from zodi_model.get_zmod import get_zmod
import numpy as np

# Single observation
result = get_zmod(
    lambda_=1.25,           # wavelength (microns)
    phase_type='kelsall',   # 'kelsall' or 'skysurf'
    day=100.0,              # day since 1990.0
    lon=90.0,               # ecliptic longitude (deg)
    lat=0.0                 # ecliptic latitude (deg)
)
print(f"Intensity: {result[0]:.6f} MJy/sr")

# Array processing
days = np.linspace(1, 365, 12)
results = get_zmod(1.25, 'kelsall', days, 90.0, 0.0)
```

## API

### get_zmod(lambda_, phase_type, day, lon, lat, **kwargs)

**Parameters:**
- `lambda_` (float): Wavelength (0.2-240 μm)
- `phase_type` (str): 'kelsall' or 'skysurf'
- `day` (float/array): Day since 1990.0
- `lon` (float/array): Ecliptic longitude (degrees)
- `lat` (float/array): Ecliptic latitude (degrees)
- `zpar` (array, optional): Custom parameters (256 elements)
- `solar_irr` (float, optional): Solar irradiance (in MJy/sr) corresponding to bandpass
- `new_iso_comp` (bool, optional): If True, include isotropic component (default: False)
- `iso_comp_only` (bool, optional): If True, model ONLY isotropic component (default: False)

**Returns:** `numpy.ndarray` - Intensity in MJy/sr

### read_zpars(filename=None)
Load model parameters from JSON.

### solar_sp(wavelength)
Solar spectrum estimate (0.2-4.85 μm).

## Model Components

1. **Smooth Cloud** - Main interplanetary dust
2. **Circumsolar Ring** - Dense dust near Sun
3. **Dust Bands** - Four asteroid bands
4. **Resonant Ring** - Earth's 1:1 orbital resonance
5. **Blobs** - Leading/trailing features

## Validation

- **Kelsall**: 0.0002% error vs IDL
- **SKYSURF**: 0.0005% error vs IDL

Test: λ=1.25 μm, day=1.0, lon=0°, lat=90° → 0.114937 MJy/sr

## Citation

```bibtex

@ARTICLE{OBrien2025,
   author = {{O'Brien}, R. and {Arendt}, R.~G. and {Windhorst}, R.~A. and {Acharya}, T. and others},
   title = "{SKYSURF-11: A New Zodiacal Light Model Optimized for Optical Wavelengths}",
   journal = {arXiv, astro-ph},
   year = 2025,
   eprint = {2510.18231},
   doi = {10.48550/arXiv.2510.18231}
}
```

## License

MIT License

## Authors

Rosalia O'Brien, Tejovrash Acharya

---

For detailed examples, see `ZodiModel_Tutorial.ipynb`
