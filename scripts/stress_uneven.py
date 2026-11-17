"""Stress for configurations with uneven parameter spacings.

core2.Problem hard-codes the uniform gaps (loggap = log|i-j|).  Only loggap has to change:
with parameters x, loggap_ij = log(x_j - x_i).  Everything else -- grad_rows, the SLP, the
multiplier LP -- then applies unchanged, so an uneven configuration can be polished to a KKT
point and its stress computed exactly as for the uniform ladder.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from core2 import Problem            # noqa: E402
from slp2 import slp, stress          # noqa: E402
import align_core as core             # noqa: E402


def problem_for(x, gauge="v2"):
    m = len(x)
    P = Problem(m, gauge=gauge)
    x = np.asarray(x, float)
    P.loggap = np.log(x[P.ju] - x[P.iu])
    return P


def to_v2(v, x):
    """regauge to v0=0, v1=1 (the Problem's gauge); the parameters are untouched"""
    v = np.asarray(v, complex)
    return (v - v[0]) / (v[1] - v[0]), np.asarray(x, float)


def load(name):
    z = np.load(ROOT / name)
    return z["v"], z["x"]


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="k4")
    ap.add_argument("--ms", default="200,400")
    args = ap.parse_args()
    for m in [int(s) for s in args.ms.split(",")]:
        path = (f"results/uneven_k4_m{m}.npz" if args.tag == "k4"
                else f"results/uneven_randhom_0_m{m}.npz")
        v, x = load(path)
        v, x = to_v2(v, x)
        P = problem_for(x)
        print(f"\n{args.tag} m={m}: D = {np.exp(P.Lvec(v).max() - P.Lvec(v).min()):.9f}")
        st = stress(P, v, tol=1e-8)
        print(f"  as stored: KKT status {st['status']} "
              f"(tight {len(st['pairs_max'])}+{len(st['pairs_min'])})", flush=True)
        v2, r = slp(P, v, rho0=1e-3, iters=60)
        st2 = stress(P, v2, tol=1e-8)
        print(f"  after uneven SLP: D = {np.exp(P.Lvec(v2).max() - P.Lvec(v2).min()):.9f}, "
              f"status {st2['status']} (tight {len(st2['pairs_max'])}+{len(st2['pairs_min'])})",
              flush=True)
        np.savez(ROOT / f"results/uneven_{args.tag}_m{m}_polished.npz", v=v2, x=x)


if __name__ == "__main__":
    main()
