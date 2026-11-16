"""Do the random-spacing optima lie on the same curve as the uniform ladder?

For each configuration, every point is projected onto the reference curve (the best known
m=800 S configuration, as a polyline in the endpoint gauge); the distance and the parameter of
the projection are recorded.  Small distances with a monotone parameter correspondence mean
"the same curve, sampled unevenly"; large distances mean a different optimum.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import align_core as core                       # noqa: E402

ALPHA = np.log(3) / np.log(4)


def gauge(v):
    w = core.gauge(np.asarray(v, complex))
    return np.conj(w) if w.imag.mean() < 0 else w


def dense(reference, N=6000):
    t = np.linspace(0, 1, N)
    return t, np.interp(t, np.linspace(0, 1, len(reference)), reference.real) + \
        1j * np.interp(t, np.linspace(0, 1, len(reference)), reference.imag)


def project(v, ref, t, R):
    """for each point of v: distance to the reference curve and the parameter there"""
    d = np.abs(v[:, None] - R[None, :])
    k = d.argmin(axis=1)
    return d[np.arange(len(v)), k], t[k]


def inversions(ti):
    return int(np.sum(np.diff(ti) < 0))


def main():
    ref = gauge(np.load(ROOT / "results/m800_converged.npy"))
    t, R = dense(ref)
    rows = json.loads((ROOT / "results/uneven-random.json").read_text())
    uni = {r["m"]: r["D"] for r in json.loads((ROOT / "results/uneven-bounded.json").read_text())
           if r["par"] == 1}
    out = []
    print(f"{'draw':>4} {'m':>4} {'D':>10} {'uniform D':>10}  {'max dist':>9} {'mean dist':>9} "
          f"{'inversions':>10}")
    for r in rows:
        z = np.load(ROOT / f"results/uneven_random_{r['par']:g}_m{r['m']}.npz")
        v = gauge(z["v"])
        d, ti = project(v, ref, t, R)
        out.append(dict(draw=int(r["par"]), m=r["m"], D=r["D"], maxdist=float(d.max()),
                        meandist=float(d.mean()), inversions=inversions(ti),
                        x=z["x"].tolist(), t=ti.tolist()))
        u = uni.get(r["m"], float("nan"))
        print(f"{int(r['par']):4d} {r['m']:4d} {r['D']:10.6f} {u:10.6f}  "
              f"{d.max():9.4f} {d.mean():9.4f} {inversions(ti):10d}")
    # control: the uniform optimum at the same m
    print("\ncontrol (uniform optima, distance to the m=800 curve):")
    for m in (16, 24, 32):
        z = np.load(ROOT / f"results/uneven_bounded_1_m{m}.npz")
        v = gauge(z["v"])
        d, ti = project(v, ref, t, R)
        out.append(dict(draw=-1, m=m, D=float(core.distortion(v)), maxdist=float(d.max()),
                        meandist=float(d.mean()), inversions=inversions(ti),
                        x=z["x"].tolist(), t=ti.tolist()))
        print(f"  uniform m={m:3d}: max dist {d.max():.4f}  mean {d.mean():.4f}  "
              f"inversions {inversions(ti)}")
    (ROOT / "results/uneven-random-analysis.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
