"""Search the m=6 problem for the least-distortion configuration.

For m=6 the problem is eight-dimensional (v_0=0, v_1=1, the other four points
complex), so the exact optimum is within reach of a wide multistart: annealed
L-BFGS on the smooth surrogate of log D, then sequential-LP polish on the exact
range, which drives the configuration to a KKT point.

    python3 scripts/optimise_m6.py --starts 3000 --polish 40

Writes data/m6_best.npy (the best configuration found, endpoint-gauged exactly
as the figure script expects) and results/m6-search.json.
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
from core2 import Problem          # noqa: E402
from slp2 import slp               # noqa: E402


def random_start(m, rng, scale=1.0):
    """Search gauge v_0 = 0, v_1 = 1 (a unit first step); the rest random.

    The chart matters: with the first step pinned to 1 the shape is anchored on a
    unit length and both the annealed L-BFGS and the LP are far better conditioned
    than in the endpoint gauge. Only the saved file is endpoint-gauged.
    """
    v = np.zeros(m, dtype=complex)
    v[1] = 1.0
    v[2:] = (rng.uniform(-0.1, 1.1, m - 2) + 1j * rng.uniform(-scale, scale, m - 2))
    return v


def anneal(P, v, betas=(1e2, 1e3, 1e4, 1e5, 1e6), maxiter=400):
    x = P.pack(v)
    for beta in betas:
        x = minimize(P.fg, x, args=(beta,), jac=True, method="L-BFGS-B",
                     options=dict(maxiter=maxiter, ftol=1e-15, gtol=1e-12)).x
    return P.unpack(x)


def gauge(v):
    return (v - v[0]) / (v[-1] - v[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--starts", type=int, default=3000)
    ap.add_argument("--polish", type=int, default=40)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--out", default="results/m6-search.json")
    args = ap.parse_args()

    m = 6
    P = Problem(m)               # search gauge v_0=0, v_1=1
    rng = np.random.default_rng(args.seed)
    t0 = time.time()

    # seed the pool with the stored configuration so a failure can never be worse
    pool = []
    stored = ROOT / "data/m6_best.npy"
    if stored.exists():
        w = P.unpack(P.pack(np.load(stored)))
        pool.append((P.D(w), w))

    collapsed = 0
    for k in range(args.starts):
        v = random_start(m, rng, scale=rng.choice([0.3, 0.6, 1.0, 1.6]))
        v = anneal(P, v)
        # the surrogate can drive two nodes together; such a start has no D to speak of
        d = np.abs(v[:, None] - v[None, :]) + np.eye(m)
        if not np.isfinite(v).all() or d.min() < 1e-3:
            collapsed += 1
            continue
        pool.append((P.D(v), v))
    print(f"  discarded {collapsed} collapsed starts", flush=True)
    pool.sort(key=lambda t: t[0])
    coarse = pool[0][0]
    print(f"  {args.starts} starts in {time.time() - t0:.0f}s: best coarse D = {coarse:.9f}", flush=True)

    # polish the distinct best ones on the exact range
    done, polished = [], []
    for D, v in pool:
        if len(polished) >= args.polish:
            break
        if any(np.abs(v - w).max() < 0.08 for w in done):
            continue
        d = np.abs(v[:, None] - v[None, :]) + np.eye(m)
        if d.min() < 1e-3:
            continue
        done.append(v)
        u, ok = v, True
        for _ in range(3):
            try:
                w, _ = slp(P, u, rho0=1e-2, iters=400)
            except ValueError:          # the LP saw a collapsed configuration
                ok = False
                break
            dw = np.abs(w[:, None] - w[None, :]) + np.eye(m)
            if not np.isfinite(w).all() or dw.min() < 1e-6:
                ok = False
                break
            u = w
        if ok:
            polished.append((P.D(u), u))
    polished.sort(key=lambda t: t[0])

    bestD, best = polished[0]
    print(f"  polished {len(polished)} distinct basins in {time.time() - t0:.0f}s")
    for D, v in polished[:8]:
        L = P.Lvec(v)
        nt = int((L >= L.max() - 1e-9).sum()), int((L <= L.min() + 1e-9).sum())
        print(f"    D = {D:.12f}   ties {nt[0]}+{nt[1]}")

    out = gauge(best)
    out = out / abs(out[-1] - out[0]) * (out[-1] - out[0]).real
    np.save(ROOT / "data/m6_best.npy", out)
    rows = dict(m=m, starts=args.starts, seed=args.seed,
                coarse=float(coarse),
                polished=[[float(d), [[float(z.real), float(z.imag)] for z in gauge(v)]]
                          for d, v in polished],
                seconds=time.time() - t0)
    p = ROOT / args.out
    p.write_text(json.dumps(rows, indent=1) + "\n")
    print(f"  wrote data/m6_best.npy (D = {bestD:.12f}) and {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
