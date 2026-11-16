"""Polish an uneven configuration at fixed parameter values, to a KKT point if possible.

Rounds of (beta-annealed L-BFGS on the smooth surrogate) followed by the exact-minimax SLP, with
the spacings frozen.  Reports D, the multiplier status and the size of the tight sets after each
round.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".ultramath/smwgx79ya/m8"))
sys.path.insert(0, str(ROOT / "scripts"))
import stress_uneven as su            # noqa: E402
import uneven_spacing as us           # noqa: E402
from slp2 import slp, stress           # noqa: E402


def to_ue(v):
    """Problem gauge (v0=0, v1=1) -> uneven_spacing gauge (v0=0, v_{m-1}=1)"""
    return v / v[-1]


def to_v2(v):
    return v / v[1]


def lbfgs_indices(P, v, x, m, betas=(1e6, 1e7, 1e8), maxiter=1000):
    """the smooth surrogate uses uneven_spacing's packing (v0=0, v_{m-1}=1), hence the conversion"""
    th = np.concatenate([to_ue(v)[1:-1].real, to_ue(v)[1:-1].imag])
    for beta in betas:
        th = minimize(lambda t: us.surrogate(np.concatenate([t, np.zeros(0)]), m, beta, 1.0,
                                             False, x),
                      th, jac=True, method="L-BFGS-B",
                      options=dict(maxiter=maxiter, ftol=1e-18, gtol=1e-14)).x
    return to_v2(us.unpack(np.concatenate([th, np.zeros(0)]), m, False, 1.0, x)[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="results/uneven_k4_m200.npz")
    ap.add_argument("--rounds", type=int, default=4)
    ap.add_argument("--rho0", type=float, default=0.05)
    ap.add_argument("--slp-iters", type=int, default=300)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    v, x = su.load(args.file)
    v, x = su.to_v2(v, x)
    m = len(x)
    P = su.problem_for(x)
    D = lambda w: float(np.exp(P.Lvec(w).max() - P.Lvec(w).min()))

    print(f"{args.file}: m={m}, frozen spacings (gap ratio {np.diff(x).max()/np.diff(x).min():.1f})")
    st = stress(P, v, tol=1e-8)
    print(f"  as stored            D = {D(v):.10f}  status {st['status']}  "
          f"tight {len(st['pairs_max'])}+{len(st['pairs_min'])}", flush=True)
    for rnd in range(args.rounds):
        v = lbfgs_indices(P, v, x, m)
        st = stress(P, v, tol=1e-8)
        print(f"  round {rnd} L-BFGS     D = {D(v):.10f}  status {st['status']}  "
              f"tight {len(st['pairs_max'])}+{len(st['pairs_min'])}", flush=True)
        v = slp(P, v, rho0=args.rho0, iters=args.slp_iters)[0]
        st = stress(P, v, tol=1e-8)
        print(f"  round {rnd} SLP        D = {D(v):.10f}  status {st['status']}  "
              f"tight {len(st['pairs_max'])}+{len(st['pairs_min'])}", flush=True)
    out = args.out or args.file.replace(".npz", "_polished.npz")
    np.savez(ROOT / out, v=v, x=x)
    print("saved", out)


if __name__ == "__main__":
    main()
