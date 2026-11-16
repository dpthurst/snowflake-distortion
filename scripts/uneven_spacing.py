"""Finite-m experiments with uneven parameter spacings.

The uniform ladder fixes the parameters x_i = i/(m-1).  Three ways to let them be uneven:

  free     no constraint at all.  Degenerate: splitting a 3-point pattern (which has D=1
           exactly, 3 points being placeable with all three ratios equal) into nearby copies
           drives D to 1 for every m, so the infimum is 1 and no curve is approximated.
  bounded  gap ratio max(dx)/min(dx) <= rho.  rho = 1 is the uniform ladder; rho > 1
           interpolates towards the degenerate one.  This is the well-posed version.
  profile  a fixed uneven profile x_i = (i/(m-1))^gamma or geometric gaps, optimising only
           the positions.  Answers "given this spacing, what shape comes out?"

    python3 scripts/uneven_spacing.py --mode bounded --rhos 1,2,4,16 --ms 8,12,16,24
    python3 scripts/uneven_spacing.py --mode profile --gammas 0.5,1,2 --ms 12,24
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
ALPHA = np.log(3) / np.log(4)


def gaps_from_w(w, m, rho):
    """m-1 positive gaps summing to 1 with max/min <= rho, and the chain-rule factor"""
    L = np.log(max(1.0 + 1e-12, rho))
    z = 0.5 * L * np.tanh(w)                      # in (-L/2, L/2)
    d = np.exp(z - z.max())
    d = d / d.sum()
    return d, 0.5 * L * (1.0 - np.tanh(w) ** 2)


def unpack(theta, m, nz, rho, x_fixed=None):
    v = np.zeros(m, dtype=complex)
    v[1:-1] = theta[:m - 2] + 1j * theta[m - 2:2 * (m - 2)]
    v[-1] = 1.0
    if x_fixed is not None:
        return v, np.asarray(x_fixed, dtype=float)
    if nz:
        w = theta[2 * (m - 2):]
        d = gaps_from_w(w, m, rho)[0]
    else:
        d = np.full(m - 1, 1.0 / (m - 1))
    x = np.concatenate([[0.0], np.cumsum(d)])
    x[-1] = 1.0
    return v, x


def exact_D(v, x):
    m = len(v)
    i, j = np.triu_indices(m, 1)
    r = np.abs(v[i] - v[j]) / (x[j] - x[i]) ** ALPHA
    return float(r.max() / r.min()), int((r >= r.max() * (1 - 1e-9)).sum()), \
        int((r <= r.min() * (1 + 1e-9)).sum())


def surrogate(theta, m, beta, rho, nz, x_fixed=None):
    v, x = unpack(theta, m, nz, rho, x_fixed)
    i, j = np.triu_indices(m, 1)
    dx = x[j] - x[i]
    d = np.abs(v[i] - v[j])
    L = np.log(d) - ALPHA * np.log(dx)
    zp, zm = beta * L, -beta * L
    mp_, mm_ = zp.max(), zm.max()
    Sp, Sm = np.exp(zp - mp_), np.exp(zm - mm_)
    val = (mp_ + np.log(Sp.sum())) / beta + (mm_ + np.log(Sm.sum())) / beta
    w = Sp / Sp.sum()
    u = Sm / Sm.sum()
    c = w - u
    G = np.zeros(m, dtype=complex)
    np.add.at(G, i, c * (v[i] - v[j]) / d ** 2)
    np.add.at(G, j, c * (v[j] - v[i]) / d ** 2)
    g = [G[1:-1].real, G[1:-1].imag]
    if nz:
        gx = np.zeros(m)
        np.add.at(gx, i, +ALPHA * c / dx)
        np.add.at(gx, j, -ALPHA * c / dx)
        dz = gaps_from_w(theta[2 * (m - 2):], m, rho)
        dg, fac = dz
        tail = np.cumsum(gx[::-1])[::-1]
        gxz = float(np.sum(gx * x))
        gz = np.array([dg[p] * (tail[p + 1] - gxz) for p in range(m - 1)])
        g.append(gz * fac)
    return val, np.concatenate(g)


def make_start(m, kind, rng):
    srcs = [ROOT / ".scratch/final_m100.npz", ROOT / "results/m200_S.npy",
            ROOT / "results/m400_converged.npy"]
    if kind < len(srcs) and srcs[kind].exists():
        z = np.load(srcs[kind])
        v = z["W"][int(np.argmin(z["Ds"]))] if hasattr(z, "files") else np.asarray(z, complex)
        t0 = np.linspace(0, 1, len(v)); t = np.linspace(0, 1, m)
        v = np.interp(t, t0, v.real) + 1j * np.interp(t, t0, v.imag)
        return (v - v[0]) / (v[-1] - v[0])
    t = np.linspace(0, 1, m)
    a, b = rng.normal(size=2)
    v = t * np.exp(1j * (a * np.sin(np.pi * t) + b * np.sin(2 * np.pi * t)))
    return (v - v[0]) / (v[-1] - v[0])


def optimise(m, rho, nz, rng, starts=5, betas=(1e2, 1e3, 1e4, 1e5, 1e6), maxiter=800):
    best = (np.inf, None, None)
    for s in range(starts):
        v0 = make_start(m, s, rng)
        th = np.concatenate([v0[1:-1].real, v0[1:-1].imag,
                             np.zeros(m - 1) if nz else np.zeros(0)])
        for beta in betas:
            res = minimize(surrogate, th, args=(m, beta, rho, nz), jac=True,
                           method="L-BFGS-B",
                           options=dict(maxiter=maxiter, ftol=1e-15, gtol=1e-12))
            th = res.x
        v, x = unpack(th, m, nz, rho)
        D, _, _ = exact_D(v, x)
        if D < best[0]:
            best = (D, v.copy(), x.copy())
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["free", "bounded", "profile", "random"], default="bounded")
    ap.add_argument("--draws", type=int, default=3)
    ap.add_argument("--rhos", default="1,2,4,16")
    ap.add_argument("--gammas", default="0.5,1,2")
    ap.add_argument("--ms", default="8,12,16,24")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    rng = np.random.default_rng(5)
    ms = [int(s) for s in args.ms.split(",")]
    rows = []
    if args.mode == "bounded":
        grid = [("bounded", float(r)) for r in args.rhos.split(",")]
    elif args.mode == "free":
        grid = [("free", np.inf)]
    elif args.mode == "random":
        grid = [("random", float(k)) for k in range(args.draws)]
    else:
        grid = [("profile", float(g)) for g in args.gammas.split(",")]
    for kind, par in grid:
        for m in ms:
            t0 = time.time()
            if kind == "profile":
                # fixed spacing, optimise positions only
                x = np.linspace(0, 1, m) ** par
                v0 = make_start(m, 0, rng)
                th = np.concatenate([v0[1:-1].real, v0[1:-1].imag])
                for beta in (1e2, 1e3, 1e4, 1e5, 1e6):
                    th = minimize(lambda t: surrogate(np.concatenate([t, np.zeros(0)]), m,
                                                      beta, 1.0, False),
                                  th, jac=True, method="L-BFGS-B",
                                  options=dict(maxiter=800, ftol=1e-15, gtol=1e-12)).x
                v, _ = unpack(np.concatenate([th, np.zeros(0)]), m, False, 1.0)
                D, nmax, nmin = exact_D(v, x)
            elif kind == "random":
                xr = np.sort(rng.uniform(size=m)); xr[0], xr[-1] = 0.0, 1.0
                v0 = make_start(m, 0, rng)
                th = np.concatenate([v0[1:-1].real, v0[1:-1].imag])
                for beta in (1e2, 1e3, 1e4, 1e5, 1e6):
                    th = minimize(lambda t: surrogate(np.concatenate([t, np.zeros(0)]), m,
                                                      beta, 1.0, False, xr),
                                  th, jac=True, method="L-BFGS-B",
                                  options=dict(maxiter=800, ftol=1e-15, gtol=1e-12)).x
                v, _ = unpack(np.concatenate([th, np.zeros(0)]), m, False, 1.0, xr)
                D, nmax, nmin = exact_D(v, xr)
                x = xr
            elif kind == "free":
                D, v, x = optimise(m, 1.0 + 1e-12, True, rng)
            else:
                D, v, x = optimise(m, par, True, rng)
            d = np.diff(x)
            rows.append(dict(mode=kind, par=par, m=m, D=D, x=x.tolist(),
                             gaps=d.tolist(), ratio=float(d.max() / d.min()),
                             seconds=time.time() - t0))
            tag = f"{kind} {par:g}"
            print(f"  {tag:14s} m={m:3d}: D={D:.9f}  gap ratio {d.max()/d.min():6.2f}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            np.savez(ROOT / f"results/uneven_{kind}_{par:g}_m{m}.npz", v=v, x=x)
    out = args.out or f"results/uneven-{args.mode}.json"
    p = Path(out) if Path(out).is_absolute() else ROOT / out
    p.write_text(json.dumps(rows, indent=1) + "\n")
    print("wrote", p)


if __name__ == "__main__":
    main()
