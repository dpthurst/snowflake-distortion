"""Repelling node sets: GUE eigenvalues, beta log-gas in a box, Riesz gas, compared with i.i.d.

Sampling is cheap for these:
  * GUE (the beta=2 log-gas): the eigenvalues of an n x n Hermitian Gaussian matrix, O(n^3) --
    exact, no MCMC.  Rescaled linearly into [0,1] they carry the semicircle density.
  * general beta / hard walls / other Riesz exponents: a plain Metropolis on the gas energy,
    O(n) per proposal, so a few million sweeps at n <= 200 is seconds.
Then the positions are optimised for each node set and D compared with the uniform and i.i.d.
families.
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
import uneven_spacing as us      # noqa: E402


def gue_nodes(n, rng):
    A = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    H = (A + A.conj().T) / np.sqrt(2)
    lam = np.linalg.eigvalsh(H)
    lam = np.sort(lam)
    return (lam - lam[0]) / (lam[-1] - lam[0])


def gas_nodes(n, rng, kernel="log", s=1.0, sweeps=4000, delta=None):
    """Metropolis for  -sum log|x_i-x_j|  (kernel='log')  or  sum |x_i-x_j|^{-s},  walls at 0,1"""
    x = np.sort(rng.uniform(size=n))
    def energy(y):
        d = np.abs(y[:, None] - y[None, :]) + np.eye(n)
        if kernel == "log":
            return -np.log(d).sum() / 2
        return (d ** (-s)).sum() / 2
    step = delta if delta else 1.0 / n
    E = energy(x)
    for _ in range(sweeps):
        moves = rng.integers(0, n, size=n)
        for i in moves:
            old = x[i]
            new = old + rng.uniform(-step, step)
            if new <= 0.0 or new >= 1.0:
                continue
            x[i] = new
            En = energy(x)
            if np.log(rng.uniform()) < (E - En):
                E = En
            else:
                x[i] = old
        acc_probe = None
    return np.sort(x)


def optimise_positions(m, x, v0, betas=(1e2, 1e3, 1e4, 1e5, 1e6, 1e7), maxiter=700):
    th = np.concatenate([v0[1:-1].real, v0[1:-1].imag])
    for beta in betas:
        th = minimize(lambda t: us.surrogate(np.concatenate([t, np.zeros(0)]), m, beta, 1.0,
                                             False, x),
                      th, jac=True, method="L-BFGS-B",
                      options=dict(maxiter=maxiter, ftol=1e-16, gtol=1e-13)).x
    v, _ = us.unpack(np.concatenate([th, np.zeros(0)]), m, False, 1.0, x)
    return v


def start_from_ref(m, ref):
    t = np.linspace(0, 1, m)
    v = np.interp(t, np.linspace(0, 1, len(ref)), ref.real) + \
        1j * np.interp(t, np.linspace(0, 1, len(ref)), ref.imag)
    return (v - v[0]) / (v[-1] - v[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ms", default="32,64")
    ap.add_argument("--draws", type=int, default=2)
    ap.add_argument("--out", default=str(ROOT / "results/repelling-nodes.json"))
    args = ap.parse_args()
    rng = np.random.default_rng(4)
    ref = np.load(ROOT / "results/m800_converged.npy")
    import align_core as core
    ref = core.gauge(ref)
    rows = []
    for m in [int(s) for s in args.ms.split(",")]:
        row = {}
        xu = np.linspace(0, 1, m)
        vu = optimise_positions(m, xu, start_from_ref(m, ref))
        row["D_uniform"] = us.exact_D(vu, xu)[0]
        for name, gen in (("gue", lambda: gue_nodes(m, rng)),
                          ("loggas", lambda: gas_nodes(m, rng, "log", sweeps=3000)),
                          ("riesz1", lambda: gas_nodes(m, rng, "riesz", s=1.0, sweeps=3000)),
                          ("riesz2", lambda: gas_nodes(m, rng, "riesz", s=2.0, sweeps=3000)),
                          ("iid", lambda: np.sort(rng.uniform(size=m)))):
            Ds, ratios = [], []
            for k in range(args.draws):
                x = gen(); x[0], x[-1] = 0.0, 1.0
                v = optimise_positions(m, x, start_from_ref(m, ref))
                Ds.append(us.exact_D(v, x)[0])
                d = np.diff(x); ratios.append(d.max() / d.min())
            row[name] = Ds
            row[name + "_ratio"] = ratios
            print(f"  m={m:4d} {name:8s}: D = " + ", ".join(f"{d:.6f}" for d in Ds) +
                  f"   gap ratio {min(ratios):.1f}-{max(ratios):.1f}", flush=True)
        print(f"  m={m:4d} uniform : D = {row['D_uniform']:.6f}", flush=True)
        rows.append(row)
    Path(args.out).write_text(json.dumps(rows, indent=1) + "\n")


if __name__ == "__main__":
    main()
