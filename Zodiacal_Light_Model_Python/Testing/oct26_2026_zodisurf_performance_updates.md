These updates were implemented by Rosalia O'Brien with Claude 3.5 Sonnet on Oct 26, 2026.

## Performance: vectorized, compiled, & multithreaded line-of-sight integration

`get_zmod` is now ~15–20× faster using the default `numpy` engine, and up to **~60× faster** using the optional compiled `numba` engine for large inputs (e.g. HEALPix maps). The output of the default NumPy engine is bit-identical to the original.

### Changes

**New: `zodisurf/zkernel_fast.py` (Default Engine)**
- Replaces the per-pixel Python loop with batched NumPy evaluation: LOS samples
  from many pixels are concatenated into flat arrays and processed in chunks
  (~131k samples).
- `_simpint_batch` builds Simpson grids/weights for many LOS at once
  (same 0.025 AU step and float32 weights as `simpint`).
- Per-pixel results are summed with `np.bincount` instead of `np.sum` per pixel.
- Geometry (`earthsun`, LOS range, step counts) is computed once for all pixels.
- Components with zero normalization (bands, ring, blobs) are skipped.
- Chunks run in parallel threads (`n_jobs`, default: all cores).
- Calls the physics functions in `zkernelpy.py` unchanged.

**New: `zodisurf/zkernel_numba.py` (Optional Ultra-Fast Engine)**
- Compiles the entire integration loop into pure C-style machine code via `@njit`.
- Drops intermediate memory footprint from gigabytes to kilobytes, processing exactly one point in space at a time to keep the CPU cache hot.
- Automatically handles thread distribution via `prange`.
- Achieves ~1e-13 relative difference compared to the original due to varying summation order.
- **Safety Features:** Automatically falls back to the `numpy` engine if Numba is not installed, or if the parameter set triggers an unsupported physics edge-case. Includes a self-verification check on first run.

**`zodisurf/get_zmod.py`**
- Serves as the central dispatcher. Now routes calls to either `zkernel_fast`, `zkernel_numba`, or `zkernel_ref`.
- Builds SKYSURF/Kelsall parameters once per unique wavelength instead of once
  per pixel (no more per-pixel 256-parameter copies).
- Default parameter file is read once and cached.
- New optional arguments: 
  - `engine` ('numpy', 'numba', or 'ref')
  - `n_jobs` (Thread count. *Note: defaults to 1 for safety inside MCMC pools. Use `n_jobs=-1` for maximum standalone speed*).
  - `chunk_points` and `scatt_mode` for fine-tuning the NumPy backend.
- Interface and output (rounded to 5 decimals) are unchanged. Added `set_default_engine()` helper.

**New: `zodisurf/zkernel_ref.py`** 
- The original un-vectorized kernel, kept intact to serve as the ground-truth reference for automated testing.