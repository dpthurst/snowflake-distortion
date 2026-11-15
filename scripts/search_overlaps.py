"""Search one configuration for pairs of index intervals with similar shape.

For intervals I=[a,a+L] and J=[c,c+K] we ask whether
    V(c + K t)  ~  s (V(a + L t)) + d,   t in [0,1],
with s, d complex — a similarity carrying the block onto the arc, with the
parameter taken affinely.  This is exactly the page's copy-mode model, so every
hit can be opened there.

The score is scale free: both sampled curves are centred and divided by their
own diameter, and only a rotation is fitted.  A small residual therefore means
"the same shape", independent of the size ratio, which is reported separately.
(Scoring the raw least-squares fit instead would let a small arc match almost
anything, because its errors are never amplified.)

Four parameters (a, L, c, K) rule out a full search, so a coarse grid over
L, K/L, a and c with 24 sample points is followed by a local refit of all four
at 201 points.

    python3 scripts/search_overlaps.py [curve.npy] [--out overlaps.json]
    python3 scripts/search_overlaps.py --surrogate 3      # random-turn control
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import align_core as core  # noqa: E402

DEFAULT = ROOT / ".ultramath/smwhdqmht/m1/s_m1600_plain.npy"


def sample(V, pos):
    pos = np.clip(pos, 0, len(V) - 1)
    return np.interp(pos, np.arange(len(V)), V.real) + \
           1j * np.interp(pos, np.arange(len(V)), V.imag)


def shape_score(block, arc):
    """Scale-free shape distance; returns (rotation, residual, size ratio).

    The rotation s satisfies arc/diam(arc) ~ s * block/diam(block), so
    size*s is the similarity carrying the block onto the arc.
    """
    bc, ac = block - block.mean(), arc - arc.mean()
    dw = float(np.abs(bc).max()) or 1.0
    da = float(np.abs(ac).max()) or 1.0
    wn, an = bc / dw, ac / da
    s = np.vdot(wn, an) / max(np.vdot(wn, wn).real, 1e-300)
    return s, float(np.abs(an - s * wn).max()), da / dw


def coarse(V, lengths, ratios, ts, step_a, step_c, thresh):
    m = len(V)
    hits = []
    for L in lengths:
        for ratio in ratios:
            K = int(round(ratio * L))
            if K < 20 or K > m - 1 or L > m - 1:
                continue
            cgrid = np.arange(0, m - K, step_c)
            arcs = sample(V, cgrid[:, None] + ts[None, :] * K)
            ac0 = arcs - arcs.mean(axis=1, keepdims=True)
            da = np.abs(ac0).max(axis=1)
            ac = ac0 / da[:, None]
            den = (ac.conj() * ac).real.sum(axis=1)
            for a in range(0, m - L, step_a):
                w = sample(V, a + ts * L)
                wc0 = w - w.mean()
                dw = float(np.abs(wc0).max())
                if dw <= 0:
                    continue
                wc = wc0 / dw
                s = (ac.conj() * wc[None, :]).sum(axis=1) / den
                resid = np.abs(wc[None, :] - s[:, None] * ac).max(axis=1)
                for k in np.flatnonzero(resid < thresh):
                    size = float(da[k] / dw)
                    if abs(size - 1) < 0.1 and abs(np.degrees(np.angle(s[k]))) < 6:
                        continue                     # too close to the identity
                    hits.append((float(resid[k]), int(a), int(L), int(cgrid[k]), int(K)))
    return hits


def refine(V, hit, nt=201, iters=70):
    """local refit of a, L, c, K by coordinate descent on the shape residual"""
    _, a, L, c, K = hit
    x = np.array([float(a), float(L), float(c), float(K)])
    t = np.linspace(0, 1, nt)
    n = len(V) - 1

    def cost(y):
        a_, L_, c_, K_ = y
        if L_ < 20 or K_ < 20 or a_ < 0 or c_ < 0 or a_ + L_ > n or c_ + K_ > n:
            return (9e9, 0j, 0.0)
        w = sample(V, a_ + t * L_)
        q = sample(V, c_ + t * K_)
        s_, r_, size_ = shape_score(w, q)      # (rotation, residual, size ratio)
        return r_, s_, size_

    val, s, size = cost(x)
    step = np.array([8.0, 8.0, 8.0, 8.0])
    for _ in range(iters):
        moved = False
        for i in range(4):
            for sgn in (+1, -1):
                y = x.copy()
                y[i] += sgn * step[i]
                v2, s2, z2 = cost(y)
                if v2 < val:
                    x, val, s, size = y, v2, s2, z2
                    moved = True
        if not moved:
            step *= 0.5
            if step.max() < 1e-3:
                break
    return dict(residual=val, a=float(x[0]), L=float(x[1]), c=float(x[2]), K=float(x[3]),
                size=float(size), argC=float(np.degrees(np.angle(s))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("curve", nargs="?", default=str(DEFAULT))
    ap.add_argument("--out", default=None)
    ap.add_argument("--step-a", type=int, default=12)
    ap.add_argument("--step-c", type=int, default=12)
    ap.add_argument("--thresh", type=float, default=0.12)
    ap.add_argument("--surrogate", type=int, default=0,
                    help="random-turn walk with the same step lengths")
    ap.add_argument("--shuffle-turns", type=int, default=0,
                    help="same step lengths and the same multiset of turns, order shuffled")
    args = ap.parse_args()

    V = core.gauge(np.load(args.curve))
    if args.shuffle_turns:
        rng = np.random.default_rng(args.shuffle_turns)
        steps = np.diff(V)
        turns = np.angle(steps[1:] / steps[:-1])
        rng.shuffle(turns)
        ph = np.concatenate([[np.angle(steps[0])], np.cumsum(turns) + np.angle(steps[0])])
        V = core.gauge(np.concatenate([[0], np.cumsum(np.abs(steps) * np.exp(1j * ph))]))
    if args.surrogate:
        rng = np.random.default_rng(args.surrogate)
        steps = np.abs(np.diff(V))
        ph = rng.uniform(0, 2 * np.pi) + np.cumsum(rng.uniform(-np.pi, np.pi, len(steps)))
        V = core.gauge(np.concatenate([[0], np.cumsum(steps * np.exp(1j * ph))]))
    m = len(V)
    ts = np.linspace(0, 1, 24)
    lengths = [120, 180, 240, 320, 400, 520, 640]
    ratios = [0.5, 0.7, 0.85, 1.0, 1.2, 1.45, 1.8, 2.2]
    tag = (f"surrogate {args.surrogate}" if args.surrogate else
           f"shuffled turns {args.shuffle_turns}" if args.shuffle_turns else "as stored")
    print(f"{Path(args.curve).name} [{tag}]: {m} nodes, coarse {len(lengths)}x{len(ratios)}, "
          f"steps {args.step_a}/{args.step_c}", flush=True)
    hits = coarse(V, lengths, ratios, ts, args.step_a, args.step_c, args.thresh)
    print(f"coarse hits below {args.thresh}: {len(hits)}", flush=True)

    hits.sort()
    chosen = []
    for h in hits:
        if all(abs(h[1] - g[1]) > 20 or abs(h[3] - g[3]) > 20 or abs(h[2] - g[2]) > 20 for g in chosen):
            chosen.append(h)
        if len(chosen) >= 150:
            break
    out = []
    for h in chosen:
        r = refine(V, h)
        a0, L, c0, K = r["a"], r["L"], r["c"], r["K"]
        dev = max(abs(c0 + K * t - (a0 + L * t)) for t in (0.0, 1.0))
        if dev < 0.1 * L and abs(r["size"] - 1) < 0.1 and abs(r["argC"]) < 4:
            continue                                  # the identity on the overlap
        r["map_deviation"] = float(dev)
        r["scale_over_alpha"] = float(r["size"] / (K / L) ** core.ALPHA)
        out.append(r)
    out.sort(key=lambda r: r["residual"])
    print(f"\n{'resid':>7s} {'block':>15s} {'arc':>15s} {'size':>7s} {'arg':>8s} {'ov':>5s} "
          f"{'size/rho^a':>10s}")
    for r in out[:20]:
        a0, a1 = r["a"], r["a"] + r["L"]
        c0, c1 = r["c"], r["c"] + r["K"]
        ov = max(0.0, min(a1, c1) - max(a0, c0))
        print(f"{r['residual']:7.4f} [{a0:6.0f},{a1:6.0f}] [{c0:6.0f},{c1:6.0f}] "
              f"{r['size']:7.3f} {r['argC']:+8.1f} {ov:5.0f} {r['scale_over_alpha']:10.3f}")
    for r in out:
        for k, v in list(r.items()):
            if isinstance(v, complex):
                print("complex field:", k, v)
                r[k] = [v.real, v.imag]
            elif isinstance(v, np.floating):
                r[k] = float(v)
            elif isinstance(v, np.integer):
                r[k] = int(v)
    if args.out:
        Path(args.out).write_text(json.dumps(out[:60], indent=1) + "\n")
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
