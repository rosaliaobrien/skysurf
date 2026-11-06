#!/usr/bin/env python3
"""
Apply Tim's and Mischa's diffs to separate copies of the zodisurf repo, run
get_zmod in a fresh interpreter for each copy across several test cases, and
check that the outputs agree.

Usage:
    python test_diff_agreement.py --repo /path/to/repo_root \
        --tim tim_zkernel.diff --mischa mischa_zkernel.diff
    # add --rtol 1e-12 if the comparison should allow floating-point noise

--repo is the directory that CONTAINS the zodisurf/ package.
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


# ----------------------------------------------------------------------
# Patching
# ----------------------------------------------------------------------
def make_patched_copy(repo: Path, diff: Path, workdir: Path, label: str, strip: int) -> Path:
    """Copy the whole repo and apply `diff` to the copy. Returns the copy's root."""
    dest = workdir / label
    shutil.copytree(repo, dest, ignore=shutil.ignore_patterns("__pycache__", ".git"))

    cmd = ["patch", f"-p{strip}", "-d", str(dest), "-i", str(diff.resolve())]

    dry = subprocess.run(cmd + ["--dry-run"], capture_output=True, text=True)
    if dry.returncode != 0:
        sys.exit(f"[{label}] dry run failed for {diff}:\n{dry.stdout}{dry.stderr}")

    real = subprocess.run(cmd, capture_output=True, text=True)
    report = real.stdout + real.stderr
    if real.returncode != 0:
        sys.exit(f"[{label}] patch failed for {diff}:\n{report}")

    # Fuzz or offsets mean the diff did not match the original exactly.
    for line in report.splitlines():
        if "fuzz" in line or "offset" in line:
            print(f"[{label}] WARNING: {line}")
    return dest


# ----------------------------------------------------------------------
# Worker: runs in its own interpreter, against one patched copy
# ----------------------------------------------------------------------
def build_cases(n):
    """Test cases for get_zmod. Phase types and fitting options cover the main code paths."""
    fit = dict(
        fitting_hong=np.array([0.5, 0.3, 0.2, 0.4, 0.5, 0.1]),
        fitting_albedo=np.linspace(0.2, 0.4, n),
        fitting_em_smooth=0.05,
        fitting_em_band=np.linspace(0.01, 0.03, n),
        fitting_em_ring=0.02,
    )
    solar = np.full(n, 200.0)
    return {
        "kelsall": dict(phase_type="kelsall"),
        "skysurf": dict(phase_type="skysurf"),
        "skysurf_iso": dict(phase_type="skysurf", new_iso_comp=True),
        "skysurf_iso_only": dict(phase_type="skysurf", new_iso_comp=True, iso_comp_only=True),
        "fitting": dict(phase_type="fitting", solar_irr=solar, **fit),
        "fitting_L2_iso": dict(phase_type="fitting", solar_irr=solar, L2=True,
                               new_iso_comp=True, **fit),
    }


def run_worker(pkg_root: Path, out_file: Path, engines):
    pkg_root = Path(pkg_root).resolve()

    # Put the zodisurf/ directory on the path so `import get_zmod` finds this copy's file.
    # get_zmod.py uses absolute `zodisurf.` imports, so the repo root must be on the path too.
    sys.path.insert(0, str(pkg_root / "zodisurf"))
    sys.path.insert(0, str(pkg_root))

    from get_zmod_old import get_zmod as gz

    # Make sure we loaded the patched copy and not an installed or cached version.
    loaded = Path(gz.__file__).resolve()
    if not loaded.is_relative_to(pkg_root):
        sys.exit(f"Loaded get_zmod from {loaded}, expected a file under {pkg_root}")

    get_zmod = gz.get_zmod

    lam = np.array([0.5, 1.25, 2.2, 3.0])
    day = np.array([3650.0, 3700.5, 4000.0, 4200.25])
    lon = np.array([10.0, 95.0, 180.0, 300.0])
    lat = np.array([0.0, 15.0, -30.0, 60.0])

    results = {}
    for engine in engines:
        for name, kw in build_cases(lam.size).items():
            key = f"{engine}/{name}"
            kw = dict(kw)
            phase_type = kw.pop("phase_type")
            try:
                out = get_zmod(lam, phase_type, day, lon, lat, engine=engine, **kw)
                results[key] = np.asarray(out, dtype=float)
            except Exception as e:
                results[key] = np.array(f"ERROR {type(e).__name__}: {e}")
    np.savez(out_file, **results)


# ----------------------------------------------------------------------
# Comparison
# ----------------------------------------------------------------------
def compare(a_file: Path, b_file: Path, rtol: float, atol: float) -> bool:
    a = np.load(a_file, allow_pickle=False)
    b = np.load(b_file, allow_pickle=False)
    ok_all = True

    if set(a.files) != set(b.files):
        print("FAIL: different case sets:", set(a.files) ^ set(b.files))
        return False

    for key in sorted(a.files):
        x, y = a[key], b[key]
        if x.dtype.kind in "US" or y.dtype.kind in "US":
            ok = x.dtype.kind == y.dtype.kind and np.array_equal(x, y)
            print(f"{'PASS' if ok else 'FAIL'} {key:30s} "
                  f"{'same error' if ok else f'{x} vs {y}'}")
        elif x.shape != y.shape:
            ok = False
            print(f"FAIL {key:30s} shape {x.shape} vs {y.shape}")
        elif np.array_equal(x, y, equal_nan=True):
            ok = True
            print(f"PASS {key:30s} bit-identical")
        else:
            ok = np.allclose(x, y, rtol=rtol, atol=atol, equal_nan=True)
            diff = np.nanmax(np.abs(x - y))
            print(f"{'PASS' if ok else 'FAIL'} {key:30s} max|diff|={diff:.3e} "
                  f"(rtol={rtol}, atol={atol})")
        ok_all &= ok
    return ok_all


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, help="Directory containing the zodisurf/ package")
    ap.add_argument("--tim", type=Path, help="Tim's diff")
    ap.add_argument("--mischa", type=Path, help="Mischa's diff")
    ap.add_argument("--strip", type=int, default=1, help="patch -p level (default 1)")
    ap.add_argument("--engines", nargs="+", default=["numpy"], choices=["numpy", "numba"])
    ap.add_argument("--rtol", type=float, default=0.0,
                    help="0 = bit-identical required; try 1e-12 for numba vs numpy")
    ap.add_argument("--atol", type=float, default=0.0)
    ap.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--pkg-root", type=Path, help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.worker:
        run_worker(args.pkg_root, args.out, args.engines)
        return 0

    if not (args.repo and args.tim and args.mischa):
        ap.error("--repo, --tim and --mischa are required")

    with tempfile.TemporaryDirectory() as td:
        workdir = Path(td)
        roots = {
            "tim": make_patched_copy(args.repo, args.tim, workdir, "tim", args.strip),
            "mischa": make_patched_copy(args.repo, args.mischa, workdir, "mischa", args.strip),
        }

        outs = {}
        for label, root in roots.items():
            outs[label] = workdir / f"{label}_results.npz"
            cmd = [sys.executable, str(Path(__file__).resolve()), "--worker",
                   "--pkg-root", str(root), "--out", str(outs[label]),
                   "--engines", *args.engines]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                sys.exit(f"[{label}] worker crashed:\n{r.stdout}{r.stderr}")

        print(f"Comparing Tim vs Mischa over engines={args.engines}\n")
        ok = compare(outs["tim"], outs["mischa"], args.rtol, args.atol)

    print("\nRESULT:", "PASS: diffs agree" if ok else "FAIL: diffs disagree")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())