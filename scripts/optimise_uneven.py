"""The finite-m problem with free parameter spacings.

The ladder of the note samples the arc at the m uniform parameters.  Here the parameters
0=x_0<...<x_{m-1}=1 are free as well:

    D = min over v in C^m and x in the simplex of  max_{i<j} r_ij / min_{i<j} r_ij,
        r_ij = |v_i - v_j| / (x_j - x_i)^alpha.

Gauge: v_0=0, v_{m-1}=1 (similarity) and x_0=0, x_{m-1}=1 (affine reparametrisation).
Unknowns: v_1..v_{m-2} complex and the m-2 free log-increments z (softmax over the m-1 gaps).

Optimised by annealed L-BFGS on the smooth surrogate of log D, with an analytic gradient.

    python3 scripts/optimise_uneven.py --ms 6,8,10,12,16,20 [--uniform]
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


def unpack(theta, m):
    """theta -> (v (m,), x (m,)) with the gauges v0=0, v_{m-1}=1, x0=0, x_{m-1}=1"""
    v = np.zeros(m, dtype=complex)
    v[1:-1] = theta[:m - 2] + 1j * theta[m - 2:2 * (m - 2)]
    v[-1] = 1.0
    z = np.concatenate([theta[2 * (m - 2):], [0.0]])       # m-1 log-gaps, last fixed at 0
    e = np.exp(z - z.max())
    d = e / e.sum()
    x = np.concatenate([[0.0], np.cumsum(d)])[:m]
    x[-1] = 1.0
    return v, x


def exact_D(v, x):
    m = len(v)
    i, j = np.triu_indices(m, 1)
    r = np.abs(v[i] - v[j]) / (x[j] - x[i]) ** ALPHA
    return float(r.max() / r.min()), int((r >= r.max() * (1 - 1e-9)).sum()), \
        int((r <= r.min() * (1 + 1e-9)).sum())


def residuals(v, x):
    m = len(v)
    i, j = np.triu_indices(m, 1)
    dx = x[j] - x[i]
    d = np.abs(v[i] - v[j])
    L = np.log(d) - ALPHA * np.log(dx)
    return i, j, d, dx, L


def surrogate(theta, m, beta, grad=True):
    v, x = unpack(theta, m)
    i, j, d, dx, L = residuals(v, x)
    zmax = beta * L
    zmin = -beta * L
    mx = zmax.max(); mn = zmin.max()
    Sp = np.exp(zmax - mx); Sm = np.exp(zmin - mn)
    val = (mx + np.log(Sp.sum())) / beta + (mn + np.log(Sm.sum())) / beta
    if not grad:
        return val
    w = Sp / Sp.sum(); u = Sm / Sm.sum()
    c = w - u                                    # dF/dL_ij up to the 1/beta scaling
    # d/d vbar_i of |v_i-v_j|^2 is (v_i-v_j): the real gradient is 2*c*(v_i-v_j)/d^2 summed
    G = np.zeros(m, dtype=complex)
    np.add.at(G, i, c * (v[i] - v[j]) / d ** 2)
    np.add.at(G, j, c * (v[j] - v[i]) / d ** 2)
    gv = np.concatenate([G[1:-1].real, G[1:-1].imag])
    # d/dx_k: L_ij has -alpha log(x_j - x_i)
    gx = np.zeros(m)
    np.add.at(gx, i, +ALPHA * c / dx)          # dL_ij/dx_i = +alpha/(x_j-x_i)
    np.add.at(gx, j, -ALPHA * c / dx)          # dL_ij/dx_j = -alpha/(x_j-x_i)
    # chain rule through the softmax of the m-1 gaps: x_k = sum_{p<k} d_p,
    # d_p = e^{z_p}/S with z_{m-2} = 0 fixed, so dx_k/dz_p = d_p (1[p<k] - x_k).
    d = np.diff(x)                             # d_0..d_{m-2}
    e = d / d.sum()
    tail = np.cumsum(gx[::-1])[::-1]           # tail[p] = sum_{k>=p} gx_k
    gxz = float(np.sum(gx * x))
    gz = np.array([e[p] * (tail[p + 1] - gxz) for p in range(m - 2)])
    g = np.concatenate([gv, gz])
    return val, g


def optimise(m, starts, betas=(1e2, 1e3, 1e4, 1e5, 1e6), maxiter=800, uniform=False,
             verbose=False):
    best = (np.inf, None, None)
    for k, (v0, x0) in enumerate(starts):
        z = np.zeros(m - 1)
        if not uniform:
            d = np.diff(x0)
            z[:-1] = np.log(np.maximum(d[:-1] / d[-1], 1e-6))
        theta = np.concatenate([v0[1:-1].real, v0[1:-1].imag,
                                z[:-1] if not uniform else np.zeros(m - 2)])
        free = theta.copy()
        for beta in betas:
            if uniform:
                n = 2 * (m - 2)
                f = lambda t: surrogate(np.concatenate([t, np.zeros(m - 2)]), m, beta)
                res = minimize(f, free[:n], jac=True, method="L-BFGS-B",
                               options=dict(maxiter=maxiter, ftol=1e-15, gtol=1e-12))
                free = np.concatenate([res.x, np.zeros(m - 2)])
            else:
                res = minimize(surrogate, free, args=(m, beta), jac=True, method="L-BFGS-B",
                               options=dict(maxiter=maxiter, ftol=1e-15, gtol=1e-12))
                free = res.x
            v, x = unpack(free, m)
            if verbose:
                print(f"    start {k} beta {beta:8.0f}: {exact_D(v, x)[0]:.9f}", flush=True)
        v, x = unpack(free, m)
        D, nmax, nmin = exact_D(v, x)
        if D < best[0]:
            best = (D, v.copy(), x.copy())
    return best


def resample(v, m):
    t0 = np.linspace(0, 1, len(v))
    t = np.linspace(0, 1, m)
    return np.interp(t, t0, v.real) + 1j * np.interp(t, t0, v.imag)


def starts_for(m, rng):
    out = []
    srcs = [ROOT / ".scratch/final_m100.npz", ROOT / "results/m200_S.npy",
            ROOT / "results/m400_converged.npy"]
    for p in srcs:
        if p.exists():
            v = np.load(p)
            if hasattr(v, "files"):
                v = v["W"][int(np.argmin(v["Ds"]))]
            out.append((resample(np.asarray(v, complex), m), np.linspace(0, 1, m)))
    out.append((np.linspace(0, 1, m).astype(complex), np.linspace(0, 1, m)))
    for _ in range(3):
        t = np.linspace(0, 1, m)
        a, b = rng.normal(size=2)
        v = t * np.exp(1j * (a * np.sin(np.pi * t) + b * np.sin(2 * np.pi * t)))
        out.append(((v - v[0]) / (v[-1] - v[0]), np.linspace(0, 1, m)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ms", default="6,8,10,12,16,20")
    ap.add_argument("--uniform", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    rng = np.random.default_rng(3)
    rows = []
    for m in [int(s) for s in args.ms.split(",")]:
        t0 = time.time()
        D, v, x = optimise(m, starts_for(m, rng), uniform=args.uniform)
        d = np.diff(x)
        rows.append(dict(m=m, D=D, gaps=d.tolist(), x=x.tolist(),
                         ratio=float(d.max() / d.min()), seconds=time.time() - t0,
                         uniform=args.uniform))
        tag = "uniform" if args.uniform else "uneven "
        print(f"  m={m:3d} {tag} D={D:.9f}  gaps in [{d.min():.2e},{d.max():.2e}] "
              f"ratio {d.max()/d.min():6.1f}  [{time.time()-t0:.0f}s]", flush=True)
        name = f"uneven_m{m}" if not args.uniform else f"uniform_m{m}"
        np.savez(ROOT / f"results/{name}.npz", v=v, x=x)
    out = Path(args.out) if args.out else ROOT / (
        "results/uneven-spacing.json" if not args.uniform else "results/uniform-spacing.json")
    if not out.is_absolute():
        out = ROOT / out
    out.write_text(json.dumps(rows, indent=1) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
