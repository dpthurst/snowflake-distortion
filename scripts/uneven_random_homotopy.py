"""i.i.d. spacings at large m, with a homotopy in the spacing.

Straight i.i.d. spacings at m >= 200 are stiff (gap ratios of 10^3-10^4), and a single L-BFGS pass
from a uniform start stalls; the spread between draws was 4%.  Here the spacings are moved
gradually from uniform to the i.i.d. target, re-optimising the positions at each step, which is
the standard continuation for this kind of family.

    python3 scripts/uneven_random_homotopy.py --ms 100,200,400 --draws 3
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


def optimise_positions(m, x, v0, betas=(1e2, 1e3, 1e4, 1e5, 1e6, 1e7), maxiter=700):
    th = np.concatenate([v0[1:-1].real, v0[1:-1].imag])
    for beta in betas:
        th = minimize(lambda t: us.surrogate(np.concatenate([t, np.zeros(0)]), m, beta, 1.0,
                                             False, x),
                      th, jac=True, method="L-BFGS-B",
                      options=dict(maxiter=maxiter, ftol=1e-16, gtol=1e-13)).x
    v, _ = us.unpack(np.concatenate([th, np.zeros(0)]), m, False, 1.0, x)
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ms", default="100,200,400")
    ap.add_argument("--draws", type=int, default=3)
    ap.add_argument("--out", default=str(ROOT / "results/uneven-random-large.json"))
    args = ap.parse_args()
    rng = np.random.default_rng(23)
    rows = []
    ref = core.gauge(np.load(ROOT / "results/m800_converged.npy"))
    for m in [int(s) for s in args.ms.split(",")]:
        for draw in range(args.draws):
            t0 = time.time()
            xt = np.sort(rng.uniform(size=m)); xt[0], xt[-1] = 0.0, 1.0
            # start from the m=800 curve sampled at the uniform grid, then walk the spacings
            t = np.linspace(0, 1, m)
            v0 = np.interp(t, np.linspace(0, 1, len(ref)), ref.real) + \
                1j * np.interp(t, np.linspace(0, 1, len(ref)), ref.imag)
            v0 = (v0 - v0[0]) / (v0[-1] - v0[0])
            xu = np.linspace(0, 1, m)
            for lam in (0.0, 0.15, 0.35, 0.6, 0.85, 1.0):
                x = (1 - lam) * xu + lam * xt
                x[0], x[-1] = 0.0, 1.0
                v0 = optimise_positions(m, x, v0)
            D = us.exact_D(v0, xt)[0]
            d = np.diff(xt)
            rows.append(dict(mode="random-homotopy", draw=draw, m=m, D=D,
                             gap_ratio=float(d.max() / d.min()),
                             seconds=time.time() - t0, x=xt.tolist()))
            np.savez(ROOT / f"results/uneven_randhom_{draw}_m{m}.npz", v=v0, x=xt)
            print(f"  draw {draw} m={m:4d}: D={D:.9f}  gap ratio {d.max()/d.min():9.1f}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            Path(args.out).write_text(json.dumps(rows, indent=1) + "\n")


if __name__ == "__main__":
    main()
