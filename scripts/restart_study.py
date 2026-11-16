"""How much does the value at fixed m depend on where the optimiser starts?

The quoted D_m is the best configuration found, hence an upper bound for the
true D_m; what is missing is an estimate of how much better a configuration we
never found could be.  This runs the same cheap restart protocol at every rung
of the ladder -- annealed L-BFGS on the smooth max/min surrogate, no SLP, no
certificate -- from a varied set of starts, and records the values reached.

At m=40 and m=100 the protocol is run with as many restarts as the earlier
census, so its cluster count and best-to-runner-up gap can be compared with the
census values (31 clusters, gap 7.3e-4 at m=40; 40 clusters, gap 7.6e-4 at
m=100).  That calibration is what licenses reading the same statistics at the
larger m, where the census cannot be afforded.

    python3 scripts/restart_study.py [--out results/restart-study.json]

Runs about half an hour in total, most of it at m=1600.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".scratch"))
sys.path.insert(0, str(ROOT / "scripts"))
from bigm import run, true_dist, V_of            # noqa: E402
import align_core as core                       # noqa: E402

ALPHA = np.log(3) / np.log(4)
BETAS = (1e3, 1e4, 1e5, 1e6)

KNOW = {                                     # best stored configuration at each m
    40: ROOT / ".scratch/final_m40.npz",
    100: ROOT / ".scratch/final_m100.npz",
    200: ROOT / "results/m200_S.npy",
    400: ROOT / "results/m400_converged.npy",
    800: ROOT / "results/m800_converged.npy",
    1600: ROOT / ".ultramath/smwhdqmht/m1/s_m1600_plain.npy",
}
ARCS = {200: ROOT / ".scratch/descended_m200.npy",
        400: ROOT / ".ultramath/smwhdqmht/m1/arc_m400.npy",
        800: ROOT / ".ultramath/smwhdqmht/m1/arc_m800.npy",
        1600: ROOT / ".ultramath/smwhdqmht/m1/arc_m1600.npy"}
RESTARTS = {40: 30, 100: 24, 200: 12, 400: 6, 800: 3, 1600: 0}


def load(path):
    if str(path).endswith(".npz"):
        z = np.load(path)
        return core.gauge(z["W"][int(np.argmin(z["Ds"]))])
    return core.gauge(np.load(path))


def gauge(v):
    v = np.asarray(v, dtype=complex)
    return (v - v[0]) / (v[1] - v[0])


def resample(v, m):
    t0 = np.linspace(0.0, 1.0, len(v))
    t = np.linspace(0.0, 1.0, m)
    return np.interp(t, t0, v.real) + 1j * np.interp(t, t0, v.imag)


def smooth_noise(m, sigma, rng, modes=4):
    """smooth complex noise vanishing at both endpoints"""
    t = np.linspace(0.0, 1.0, m)
    out = np.zeros(m, dtype=complex)
    for k in range(1, modes + 1):
        c = rng.normal(size=2)
        out += (c[0] + 1j * c[1]) * np.sin(np.pi * k * t) / k
    return sigma * out / max(np.abs(out).max(), 1e-12)


def random_curve(m, rng):
    """a random smooth open curve from 0 to 1, points spread along its length"""
    t = np.linspace(0.0, 1.0, m)
    a, b, c = rng.normal(size=3) * 0.9
    phase = a * np.sin(np.pi * t) + b * np.sin(2 * np.pi * t) + c * np.sin(3 * np.pi * t)
    return gauge(t * np.exp(1j * phase))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "results/restart-study.json"))
    ap.add_argument("--ms", default="40,100,200,400,800,1600")
    args = ap.parse_args()
    rng = np.random.default_rng(20261116)
    rows = []
    for m in [int(s) for s in args.ms.split(",")]:
        if RESTARTS[m] == 0:
            continue
        base = load(KNOW[m]) if m in KNOW and KNOW[m].exists() else None
        files = {}
        for name, path in (("known", KNOW.get(m)), ("arc", ARCS.get(m))):
            if path is not None and Path(path).exists():
                files[name] = load(path)
        print(f"\n=== m={m}, {RESTARTS[m]} restarts; starts available: "
              f"{sorted(files)}", flush=True)
        for r in range(RESTARTS[m]):
            kind = ["known+noise", "arc+noise", "coarse+noise"][r % 3]
            if kind == "known+noise":
                v0 = gauge(files["known"] + smooth_noise(m, 0.01, rng))
            elif kind == "arc+noise" and "arc" in files:
                v0 = gauge(files["arc"] + smooth_noise(m, 0.01, rng))
            elif kind == "arc+noise":                    # no arc stored: neighbour
                other = max([k for k in KNOW if k < m], default=min(KNOW))
                v0 = gauge(resample(load(KNOW[other]), m) + smooth_noise(m, 0.02, rng))
            else:                                        # resampled coarser rung
                other = max([k for k in KNOW if k < m], default=min(KNOW))
                v0 = gauge(resample(load(KNOW[other]), m) + smooth_noise(m, 0.02, rng))
            t0 = time.time()
            D, x = run(m, v0, BETAS, maxiter=400, verbose=False)
            secs = time.time() - t0
            rows.append(dict(m=m, restart=r, kind=kind, D=float(D), seconds=secs))
            print(f"  {r:2d} {kind:11s} D={D:.10f}  [{secs:6.1f}s]", flush=True)
            Path(args.out).write_text(json.dumps(rows, indent=1) + "\n")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
