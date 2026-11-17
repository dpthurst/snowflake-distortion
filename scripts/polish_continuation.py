"""Move the spacings from uniform to a target uneven set, polishing at every step.

A KKT point of the m-point minimax needs about 2m-2 tight pairs (2(m-2) equations to balance with
the two normalisations), and the certified uniform configurations have exactly that many
(48+346 at m=200, 70+706 at m=400).  An uneven configuration that has been polished *from a bad
start* does not get there: it crawls.  Moving the spacings gradually, with a polish at each step,
keeps the configuration at marginal contact throughout.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import stress_uneven as su            # noqa: E402
import uneven_spacing as us           # noqa: E402
from slp2 import slp, stress           # noqa: E402


def to_ue(v):
    return v / v[-1]


def to_v2(v):
    return v / v[1]


def polish(v, x, m, betas=(1e5, 1e6, 1e7), slp_iters=200, rho0=0.02):
    P = su.problem_for(x)
    th = np.concatenate([to_ue(v)[1:-1].real, to_ue(v)[1:-1].imag])
    for beta in betas:
        th = minimize(lambda t: us.surrogate(np.concatenate([t, np.zeros(0)]), m, beta, 1.0,
                                             False, x),
                      th, jac=True, method="L-BFGS-B",
                      options=dict(maxiter=600, ftol=1e-17, gtol=1e-14)).x
    v = to_v2(us.unpack(np.concatenate([th, np.zeros(0)]), m, False, 1.0, x)[0])
    v = slp(P, v, rho0=rho0, iters=slp_iters)[0]
    return v, P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="results/m200_S.npy")   # certified uniform m=200
    ap.add_argument("--target", default="results/uneven_k4_m200.npz")
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument("--out", default="results/uneven_k4_m200_continuation.npz")
    args = ap.parse_args()
    v = su.core.gauge(np.load(ROOT / args.start))
    m = len(v)
    v = to_v2(v)
    _, xt = su.load(args.target)
    xu = np.linspace(0, 1, m)
    print(f"m={m}: uniform certified start, target gap ratio {np.diff(xt).max()/np.diff(xt).min():.1f}")
    for lam in np.linspace(0, 1, args.steps + 1):
        x = (1 - lam) * xu + lam * xt
        x[0], x[-1] = 0.0, 1.0
        v, P = polish(v, x, m)
        st = stress(P, v, tol=1e-8)
        D = float(np.exp(P.Lvec(v).max() - P.Lvec(v).min()))
        print(f"  lambda={lam:.2f}: D = {D:.10f}  status {st['status']}  "
              f"tight {len(st['pairs_max'])}+{len(st['pairs_min'])}  "
              f"(need about {2*(m-2)})", flush=True)
        np.savez(ROOT / args.out, v=v, x=x)


if __name__ == "__main__":
    main()
