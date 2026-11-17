"""Search parameter patterns x for which the fixed-pattern m-point bound is largest.

Every curve of distortion D carries, at the parameters x, a configuration with spread at most D,
so D_inf >= D_m(x) for *every* pattern x, and hence D_inf >= sup_x D_m(x).  The uniform pattern
is only one candidate; this searches the space of patterns at small m: geometric ladders, draws
from a Dirichlet law, and coordinate-wise hill climbing from the best of those.

Two warnings, both learned the hard way.  The value returned is an *upper* bound for D_m(x) -- the
spread of an achieved configuration -- so a pattern is only useful once its D_m(x) is certified by
certify_mm.py with the same pattern.  And a failed optimisation looks exactly like an excellent
pattern: with too few starts the routine returns a value far above the true D_m(x), which a search
over patterns will happily pick as its maximum.  Hence the reported quantiles, and hence the default
of many starts.

    python3 scripts/gap_pattern_search.py --m 6 --draws 300
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
from core2 import Problem              # noqa: E402
from slp2 import slp                   # noqa: E402
from stress_uneven import problem_for  # noqa: E402


def anneal(P, v, betas=(1e2, 1e3, 1e4, 1e5, 1e6), maxiter=250):
    x = P.pack(v)
    for beta in betas:
        x = minimize(P.fg, x, args=(beta,), jac=True, method="L-BFGS-B",
                     options=dict(maxiter=maxiter, ftol=1e-15, gtol=1e-12)).x
    return P.unpack(x)


def D_of_pattern(x, rng, starts=120, polish=True, report=False):
    P = problem_for(x)
    m = len(x)
    best, bestv, vals = np.inf, None, []
    for _ in range(starts):
        v = np.zeros(m, dtype=complex)
        v[1] = 1.0
        v[2:] = rng.uniform(0, 1, m - 2) + 1j * rng.uniform(-1, 1, m - 2)
        v = anneal(P, v)
        d = np.abs(v[:, None] - v[None, :]) + np.eye(m)
        if not np.isfinite(v).all() or d.min() < 1e-9:
            continue
        if polish:
            try:
                v, _ = slp(P, v, rho0=1e-2, iters=250)
            except ValueError:
                continue
        D = P.D(v)
        if np.isfinite(D):
            vals.append(D)
            if D < best:
                best, bestv = D, v
    if report and vals:
        print(f"      starts={len(vals)}  best={best:.9f}  p5={np.percentile(vals, 5):.9f}"
              f"  median={np.percentile(vals, 50):.9f}")
    return best, bestv


def ladder(m, r):
    g = r ** np.arange(m - 1)
    x = np.concatenate([[0.0], np.cumsum(g)])
    return x / x[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=6)
    ap.add_argument("--draws", type=int, default=300)
    ap.add_argument("--starts", type=int, default=120)
    ap.add_argument("--refine", type=int, default=40)
    ap.add_argument("--min-gap", type=float, default=0.01)
    ap.add_argument("--out", default=None)
    ap.add_argument("--seed", type=int, default=2)
    args = ap.parse_args()

    m = args.m
    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    seen = {}

    def consider(x, tag):
        g = np.diff(x)
        if g.min() < args.min_gap:
            return None
        key = tuple(np.round(g / g.sum(), 4))
        if key in seen:
            return seen[key][0]
        D, v = D_of_pattern(x, rng, starts=args.starts)
        seen[key] = (D, x.copy(), tag, v)
        return D

    consider(np.linspace(0, 1, m), "uniform")
    for r in (1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0):
        consider(ladder(m, r), f"ladder r={r}")
        consider(1 - ladder(m, r)[::-1] + 0, f"ladder r={r} reversed") if False else None
    for _ in range(args.draws):
        g = rng.dirichlet(np.full(m - 1, 1.0))
        if g.min() < args.min_gap:
            continue
        consider(np.concatenate([[0.0], np.cumsum(g)]), "dirichlet")
    print(f"  {len(seen)} patterns in {time.time() - t0:.0f}s", flush=True)

    ranked = sorted(seen.values(), key=lambda t: -t[0])
    D0, x0 = ranked[0][0], ranked[0][1]
    print(f"  best so far D = {D0:.9f} gaps {np.round(np.diff(x0), 4).tolist()}  [{ranked[0][2]}]")

    # hill climb on the gaps
    best_D, best_x = D0, x0.copy()
    for it in range(args.refine):
        g = np.diff(best_x)
        scale = 0.25 * (0.97 ** it)
        trial = g * np.exp(rng.normal(0, scale, m - 1))
        trial = trial / trial.sum()
        if trial.min() < args.min_gap:
            continue
        x = np.concatenate([[0.0], np.cumsum(trial)])
        D, _ = D_of_pattern(x, rng, starts=args.starts)
        if D > best_D:
            best_D, best_x = D, x
            print(f"    refine {it}: D = {D:.9f} gaps {np.round(np.diff(x), 4).tolist()}", flush=True)
    print(f"  best pattern: D = {best_D:.9f} gaps {np.round(np.diff(best_x), 5).tolist()}")
    out = ROOT / (args.out or f"results/gap-patterns-m{m}.json")
    out.write_text(json.dumps(dict(m=m, uniform=[v for v in seen.values() if v[2] == "uniform"][0][0],
                                   best=float(best_D), x=best_x.tolist(),
                                   top=[[float(t[0]), t[1].tolist(), t[2]] for t in ranked[:20]]),
                              indent=1) + "\n")
    print(f"  wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
