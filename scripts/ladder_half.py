"""The finite-m ladder at the critical exponent alpha = 1/2.

The question flagged in the note is whether [0,1] with the distance |x-y|^{1/2} admits a
bi-Lipschitz (indeed any finite-distortion) embedding in the plane, i.e. whether
sup_m D_m(1/2) is finite.  This computes the best D_m found for m = 4 ... 64 by the same
annealed L-BFGS used for the main ladder, from a handful of starts.  It is an indicator,
not a proof in either direction.

    python3 scripts/ladder_half.py
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".scratch"))
import bigm                                     # noqa: E402

bigm.al = 0.5                                   # the critical exponent
bigm.CACHE.clear()

from bigm import run, true_dist, V_of            # noqa: E402

BETAS = (1e3, 1e4, 1e5, 1e6, 1e7)
rows = []
rng = np.random.default_rng(11)
for m in (4, 6, 8, 12, 16, 24, 32, 48, 64):
    best = np.inf
    for trial in range(6):
        t = np.linspace(0, 1, m)
        if trial == 0:                          # straight, evenly spaced
            v0 = t.astype(complex)
        elif trial == 1:                        # circular arc
            v0 = np.exp(0.9j * np.pi * (t - 0.5))
        else:                                   # smooth random
            a, b = rng.normal(size=2)
            v0 = t * np.exp(1j * (a * np.sin(np.pi * t) + b * np.sin(2 * np.pi * t)))
        v0 = (v0 - v0[0]) / (v0[1] - v0[0])
        D, x = run(m, v0, BETAS, maxiter=600, verbose=False)
        best = min(best, float(D))
    rows.append(dict(m=m, D=best))
    print(f"  m={m:3d}   best D = {best:.6f}", flush=True)
(ROOT / "results/ladder-half.json").write_text(json.dumps(rows, indent=1) + "\n")
print("wrote results/ladder-half.json")
