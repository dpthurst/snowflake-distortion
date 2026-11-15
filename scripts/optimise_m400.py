"""Converge one m=400 configuration towards a stationary point of the minimax.

The stored m=400 points are not stationary (the KKT check in
`limiting_curve_directions.md` gives relative residuals 1.07 and 1.04), and four
of them span D = 1.5423 .. 1.5526. This runs the exact-minimax SLP with a
beta-annealed L-BFGS preparation, records D and the multiplier certificate at
every round, and saves the best configuration so the overlap search can be
rerun on it.

    python3 scripts/optimise_m400.py [--m 400] [--rounds 5] [--minutes 40]

The certificate is `slp2.stress(P, v)`: it solves for multipliers on the tight
max and min pairs, so status 0 means a first-order KKT point of the minimax was
reached. Status 2 means none exists — the point is not stationary.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".scratch"))
sys.path.insert(0, str(ROOT / ".ultramath/smwgx79ya/m8"))
from bigm import phi_grad, V_of, true_dist      # noqa: E402
from core2 import Problem                       # noqa: E402
from slp2 import slp, stress                    # noqa: E402

START = ROOT / ".ultramath/smwhdqmht/m1/s_m400.npy"


def exact_D(v):
    v = np.asarray(v)
    m = len(v)
    i, j = np.triu_indices(m, 1)
    r = np.abs(v[i] - v[j]) / np.abs(i - j) ** (np.log(3) / np.log(4))
    return float(r.max() / r.min()), int((r >= r.max() * (1 - 1e-9)).sum()), \
        int((r <= r.min() * (1 + 1e-9)).sum())


def to_v2(v):
    """gauge v[0]=0, v[1]=1"""
    v = np.asarray(v, dtype=complex)
    v = (v - v[0]) / (v[1] - v[0])
    return v


def pack(v):
    return np.concatenate([v[2:].real, v[2:].imag])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=400)
    ap.add_argument("--rounds", type=int, default=5)
    ap.add_argument("--minutes", type=float, default=40.0)
    ap.add_argument("--slp-iters", type=int, default=200)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    m = args.m
    t0 = time.time()
    deadline = t0 + 60 * args.minutes
    P = Problem(m)
    v = to_v2(np.load(START))
    D, nmax, nmin = exact_D(v)
    st = stress(P, v, tol=1e-8)
    log = [dict(round=-1, D=D, ties=[nmax, nmin], stress=int(st["status"]),
                n_act=[len(st["pairs_max"]), len(st["pairs_min"])],
                seconds=0.0)]
    print(f"start: D={D:.10f} ties={nmax}+{nmin} stress={st['status']} "
          f"({len(st['pairs_max'])}+{len(st['pairs_min'])} tight pairs)", flush=True)
    best = (D, v.copy(), int(st["status"]))

    x = pack(v)
    for rnd in range(args.rounds):
        if time.time() > deadline:
            print("deadline reached", flush=True)
            break
        for beta in (1e3, 1e4, 1e5, 1e6):
            x = minimize(phi_grad, x, args=(m, beta), jac=True, method="L-BFGS-B",
                         options=dict(maxiter=150, ftol=1e-16, gtol=1e-15)).x
        v_ls = V_of(x, m)
        v_slp, _ = slp(P, v_ls, rho0=0.01, iters=args.slp_iters)
        for cand in (v_ls, v_slp):
            Dc, nmaxc, nminc = exact_D(cand)
            stc = stress(P, cand, tol=1e-8)
            log.append(dict(round=rnd, D=Dc, ties=[nmaxc, nminc],
                            stress=int(stc["status"]),
                            n_act=[len(stc["pairs_max"]), len(stc["pairs_min"])],
                            seconds=time.time() - t0))
            print(f" round {rnd}: D={Dc:.10f} ties={nmaxc}+{nminc} "
                  f"stress={stc['status']} ({len(stc['pairs_max'])}+{len(stc['pairs_min'])}) "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            if Dc < best[0] - 1e-12 or (stc["status"] == 0 and best[2] != 0):
                if Dc <= best[0] + 1e-6:
                    best = (Dc, cand.copy(), int(stc["status"]))
        x = pack(v_slp)

    D, v, status = best
    out = Path(args.out) if args.out else ROOT / f"results/m{m}_converged.npy"
    np.save(out, v)
    summary = dict(m=m, D=D, stress=status, seconds=time.time() - t0, log=log,
                   out=str(out))
    (ROOT / f"results/m{m}_converged.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(f"\nbest D={D:.10f} stress={status} -> {out}", flush=True)


if __name__ == "__main__":
    main()
