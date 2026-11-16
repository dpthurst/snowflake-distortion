"""The kappa-bounded spacing family at large m, with warm starts.

Deterministic spacings (no random draws): the gap ratio is bounded by kappa and both the
positions and the spacings are optimised.  Each m is started from the previous one
(interpolated), which is what makes m = 200, 400 affordable.

    python3 scripts/uneven_large.py --kappa 4 --ms 64,96,128,200,400
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import uneven_spacing as us                        # noqa: E402
import align_core as core                          # noqa: E402


def resample_state(v, x, m):
    t = np.linspace(0, 1, m)
    return (np.interp(t, x, v.real) + 1j * np.interp(t, x, v.imag)), \
        np.interp(t, x, np.linspace(0, 1, len(x)))


def seed_theta(v0, x0, m, kappa):
    """theta for uneven_spacing.unpack at gap ratio bound kappa, matching x0 as far as possible"""
    d = np.diff(x0)
    L = np.log(kappa)
    z = np.log(np.maximum(d / d.mean(), 1e-12))
    z = np.clip(z, -0.49 * L, 0.49 * L)
    w = np.arctanh(np.clip(2 * z / L, -0.999, 0.999)) if L > 0 else np.zeros(m - 1)
    return np.concatenate([v0[1:-1].real, v0[1:-1].imag, w])


def solve(m, kappa, theta, betas=(1e2, 1e3, 1e4, 1e5, 1e6, 1e7), maxiter=600):
    for beta in betas:
        theta = minimize(us.surrogate, theta, args=(m, beta, kappa, True), jac=True,
                         method="L-BFGS-B",
                         options=dict(maxiter=maxiter, ftol=1e-16, gtol=1e-13)).x
    v, x = us.unpack(theta, m, True, kappa)
    return theta, v, x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", type=float, default=4.0)
    ap.add_argument("--ms", default="64,96,128,200,400")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    rng = np.random.default_rng(7)
    prev = None
    rows = []
    for m in [int(s) for s in args.ms.split(",")]:
        t0 = time.time()
        # seeds: previous level, the uniform S configuration, and a random smooth curve
        seeds = []
        if prev is not None:
            v0, x0 = resample_state(prev[0], prev[1], m)
            seeds.append(seed_theta(v0, x0, m, args.kappa))
        src = ROOT / "results/m200_S.npy"
        if src.exists():
            v = core.gauge(np.load(src))
            v0, x0 = resample_state(v, np.linspace(0, 1, len(v)), m)
            seeds.append(seed_theta(v0, x0, m, args.kappa))
        v0 = us.make_start(m, 2, rng)
        seeds.append(seed_theta(v0, np.linspace(0, 1, m), m, args.kappa))
        # a deliberately uneven seed: gaps spread over the whole allowed ratio
        for _ in range(2):
            z = rng.uniform(-0.5, 0.5, size=m - 1) * np.log(args.kappa)
            d = np.exp(z - z.max()); d /= d.sum()
            xu = np.concatenate([[0.0], np.cumsum(d)])
            xu[-1] = 1.0
            vv = core.gauge(np.load(ROOT / "results/m200_S.npy"))
            v0, _ = resample_state(vv, np.linspace(0, 1, len(vv)), m)
            seeds.append(seed_theta(v0, xu, m, args.kappa))
        best = (np.inf, None, None, None)
        for th in seeds:
            th, v, x = solve(m, args.kappa, th)
            D = us.exact_D(v, x)[0]
            if D < best[0]:
                best = (D, v, x, th)
        D, v, x, th = best
        prev = (v, x)
        d = np.diff(x)
        rows.append(dict(kappa=args.kappa, m=m, D=D, seconds=time.time() - t0,
                         gaps=d.tolist(), x=x.tolist(),
                         gap_ratio=float(d.max() / d.min())))
        np.savez(ROOT / f"results/uneven_k{args.kappa:g}_m{m}.npz", v=v, x=x)
        print(f"  kappa={args.kappa:g} m={m:4d}: D={D:.9f}  gap ratio {d.max()/d.min():6.2f}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        out = Path(args.out) if args.out else ROOT / f"results/uneven-kappa{args.kappa:g}-large.json"
        if not out.is_absolute():
            out = ROOT / out
        out.write_text(json.dumps(rows, indent=1) + "\n")


if __name__ == "__main__":
    main()
