#!/usr/bin/env python
"""
compare_local.py
----------------
Compares the fast ZodiSURF model against the original un-vectorized kernel.
Run from the `Testing` directory.
"""
import sys
import os
import time
import numpy as np

# --- 1. Tell Python where to find zodisurf ---
# This finds the parent directory of this script (Zodiacal_Light_Model_Python)
# and adds it to the system path so `import zodisurf` works.
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# --- 2. Import everything ---
from zodisurf import get_zmod

# Import the original kernel directly, plus tools to prep its inputs
import zodisurf.zkernelpy_original as zold
from zodisurf.get_zmod import read_zpars, mk_zdata
from zodisurf.skysurf_params import get_albedo, get_hong_params, get_mult, get_emiss, put_zpar

ATOL = 1.5e-5
RTOL = 1e-6
DIRBE = [1.25, 2.2, 3.5, 4.9, 12., 25., 60., 100., 140., 240.]

# ======================================================================
# Reconstruct the original get_zmod behavior to feed the old kernel
# ======================================================================
def get_zmod_original(lambda_, phase_type, day, lon, lat, L2=False):
    """Mimics the original get_zmod.py to feed zkernelpy_original.py"""
    lambda_ = np.atleast_1d(np.asarray(lambda_, dtype=float))
    day = np.atleast_1d(np.asarray(day, dtype=float))
    lon = np.atleast_1d(np.asarray(lon, dtype=float))
    lat = np.atleast_1d(np.asarray(lat, dtype=float))

    if lambda_.size == 1 and lambda_.shape != day.shape:
        lambda_ = np.full_like(day, lambda_[0], dtype=float)

    zpar = read_zpars()
    
    if phase_type == 'kelsall':
        zpar = np.array([zpar]*len(lambda_), dtype=float)
        
    elif phase_type == 'skysurf':
        albedo_arr = get_albedo(lambda_)
        hong_params_arr = get_hong_params(lambda_)
        mult_arr = get_mult(lambda_)
        hong_params_arr[:,-3:] = hong_params_arr[:,3:] * mult_arr
        emiss_arr = get_emiss(lambda_)
        zpar = put_zpar(zpar, 0, 0, 0, albedo_arr, det1=0, hg3=hong_params_arr, E1=emiss_arr)

    data = mk_zdata(lambda_, day, lon, lat)
    
    # CALL THE ORIGINAL UNMODIFIED KERNEL
    zodi = zold.zkernel(data, zpar, phase_type=phase_type, L2=L2, no_colcorr=True)
    return np.round(zodi, 5)

# ======================================================================
# Setup test cases
# ======================================================================
def run_comparison():
    n = 1000
    rng = np.random.default_rng(2026)
    lon = rng.uniform(0.0, 360.0, n)
    lat = np.degrees(np.arcsin(rng.uniform(-1.0, 1.0, n)))
    day_rand = rng.uniform(1.0, 12000.0, n)
    full = lambda v, size=n: np.full(size, float(v))

    lats_scan = np.arange(-90, 90, 10).astype(float)
    
    cases = [
        ('SKYSURF 1.0um, day 200, L2',          full(1.0),  'skysurf', full(200.0), lon, lat, True),
        ('SKYSURF 1.25um, rand days, Earth',    full(1.25), 'skysurf', day_rand,    lon, lat, False),
        ('SKYSURF mixed wave, rand days, L2',   rng.choice([0.8, 1.0, 1.25, 1.6, 2.0], n), 'skysurf', day_rand, lon, lat, True),
        ('Kelsall 2.2um, day 150, Earth',       full(2.2),  'kelsall', full(150.0), lon, lat, False),
        ('Kelsall rand DIRBE, rand days, Earth',rng.choice(DIRBE, n), 'kelsall', day_rand, lon, lat, False),
        ('Notebook example: lat scan',          full(2.2, lats_scan.size), 'kelsall', full(150.0, lats_scan.size), np.zeros(lats_scan.size), lats_scan, False)
    ]

    print(f"Comparing FAST code to zkernelpy_original.py...")
    print("=" * 78)
    
    all_ok = True
    total_t_new = 0.0
    total_t_old = 0.0

    for name, lam, phase, day, ln, lt, L2 in cases:
        print(f"\n{name}   (n={lam.size})")
        
        # Run original
        t0 = time.perf_counter()
        orig_res = get_zmod_original(lam, phase, day, ln, lt, L2=L2)
        t_orig = time.perf_counter() - t0
        total_t_old += t_orig
        
        # Run new (force n_jobs=1 so threading doesn't skew small N benchmarks)
        t0 = time.perf_counter()
        new_res = get_zmod(lam, phase, day, ln, lt, L2=L2, n_jobs=1, engine='numba')
        t_new = time.perf_counter() - t0
        total_t_new += t_new
        
        print(f"  ORIGINAL: {t_orig:7.3f} s")
        print(f"  NEW     : {t_new:7.3f} s")

        # Compare
        d = np.abs(new_res - orig_res)
        bad = (d > ATOL + RTOL * np.abs(orig_res))
        ok = not np.any(bad)
        all_ok &= ok
        
        n_exact = int(np.sum(new_res == orig_res))
        print(f"  max|diff| = {np.nanmax(d):.3g}   bit-identical: {n_exact}/{lam.size}")
        print(f"  RESULT: {'PASS' if ok else 'FAIL'}")

    print("\n" + "=" * 78)
    print(f"Total compute time: ORIGINAL {total_t_old:.2f} s, NEW {total_t_new:.2f} s ({total_t_old/total_t_new:.1f}x faster)")
    print("OVERALL: " + ("PASS" if all_ok else "FAIL"))

if __name__ == '__main__':
    run_comparison()