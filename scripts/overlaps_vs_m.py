"""Partial self-overlaps as a function of m, and across m.

Three questions:

1. Does the density of partial self-overlaps change with m?  The search is run
   on a ladder with the *same relative* grid (lengths and steps as fractions of
   m), so the curves are compared on equal terms, together with shuffled-turn
   controls at each m.
2. Does an overlap found at m=800 survive at m=1600 at doubled indices?  That
   is what a limiting object would do; a finite-size artefact would not.
3. Cross-m: is a block of one curve similar to an arc of another beyond the
   trivial correspondence of the two samplings?

    python3 scripts/overlaps_vs_m.py
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import align_core as core  # noqa: E402
import search_overlaps as so  # noqa: E402

CURVES = {
    200: ROOT / "data/inputs/descended_m200.npy",
    400: ROOT / "data/inputs/s_m400.npy",
    800: ROOT / "data/inputs/s_m800_plain.npy",
    1600: ROOT / "data/inputs/s_m1600_plain.npy",
}

LENGTHS = [0.075, 0.11, 0.15, 0.20, 0.25, 0.325, 0.40]
RATIOS = [0.5, 0.7, 0.85, 1.0, 1.2, 1.45, 1.8, 2.2]


def load(m):
    return core.gauge(np.load(CURVES[m]))


def shuffled(V, seed):
    rng = np.random.default_rng(seed)
    steps = np.diff(V)
    turns = np.angle(steps[1:] / steps[:-1])
    rng.shuffle(turns)
    ph = np.concatenate([[np.angle(steps[0])], np.cumsum(turns) + np.angle(steps[0])])
    return core.gauge(np.concatenate([[0], np.cumsum(np.abs(steps) * np.exp(1j * ph))]))


def coarse2(VB, VA, mB, mA, thresh=0.2):
    """Block from VB, arc from VA; same relative grid as the within-curve search."""
    lengths = [max(24, int(round(f * (mB - 1)))) for f in LENGTHS]
    step_a = max(2, int(round((mB - 1) / 133)))
    step_c = max(2, int(round((mA - 1) / 133)))
    ts = np.linspace(0, 1, 24)
    hits = []
    for L in lengths:
        for ratio in RATIOS:
            K = int(round(ratio * L * (mA - 1) / (mB - 1)))
            if K < 20 or K > mA - 2 or L > mB - 2:
                continue
            cgrid = np.arange(0, mA - 1 - K, step_c)
            arcs = so.sample(VA, cgrid[:, None] + ts[None, :] * K)
            ac0 = arcs - arcs.mean(axis=1, keepdims=True)
            da = np.abs(ac0).max(axis=1)
            ac = ac0 / da[:, None]
            den = (ac.conj() * ac).real.sum(axis=1)
            for a in range(0, mB - 1 - L, step_a):
                w = so.sample(VB, a + ts * L)
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
                        continue
                    hits.append((float(resid[k]), int(a), int(L), int(cgrid[k]), int(K)))
    return hits


def run(VB, VA, mB, mA, thresh=0.2):
    if VB is VA:
        lengths = [max(24, int(round(f * (mB - 1)))) for f in LENGTHS]
        step = max(2, int(round((mB - 1) / 133)))
        ts = np.linspace(0, 1, 24)
        return so.coarse(VB, lengths, RATIOS, ts, step, step, thresh)
    return coarse2(VB, VA, mB, mA, thresh)


def main():
    out = {}
    print("== within-curve search, same relative grid ==")
    print(f"{'m':>5s} {'curve':>9s} {'best':>7s} {'#<0.12':>7s} {'#<0.15':>7s} {'#<0.20':>7s}")
    for m in (200, 400, 800, 1600):
        V = load(m)
        for tag, C in (("stored", V), ("shuffled", shuffled(V, 1))):
            hits = run(C, C, m, m)
            r = np.array([h[0] for h in hits])
            row = dict(best=float(r.min()) if r.size else None,
                       n12=int((r < 0.12).sum()), n15=int((r < 0.15).sum()),
                       n20=int((r < 0.20).sum()), total=int(r.size))
            out[f"within_{m}_{tag}"] = row
            print(f"{m:5d} {tag:>9s} {row['best'] if row['best'] else float('nan'):7.4f} "
                  f"{row['n12']:7d} {row['n15']:7d} {row['n20']:7d}")

    # ---- 1b. same node count: decimate every curve to 200 nodes
    print("\n== every curve decimated to 200 nodes, identical grid ==")
    print(f"{'m':>5s} {'best':>7s} {'#<0.12':>7s} {'#<0.15':>7s} {'#<0.20':>7s}")
    N = 200
    lengths = [int(round(f * (N - 1))) for f in LENGTHS]
    ts = np.linspace(0, 1, 24)
    for m in (200, 400, 800, 1600):
        V = load(m)
        idx = np.round(np.linspace(0, len(V) - 1, N)).astype(int)
        W = core.gauge(V[idx])
        hits = so.coarse(W, lengths, RATIOS, ts, 2, 2, 0.2)
        r = np.array([h[0] for h in hits]) if hits else np.array([])
        row = dict(best=float(r.min()) if r.size else None,
                   n12=int((r < 0.12).sum()), n15=int((r < 0.15).sum()),
                   n20=int((r < 0.20).sum()))
        out[f"decimated200_{m}"] = row
        print(f"{m:5d} {row['best'] if row['best'] else float('nan'):7.4f} "
              f"{row['n12']:7d} {row['n15']:7d} {row['n20']:7d}")

    # ---- 2. nesting: do m=800 overlaps reappear at m=1600?
    print("\n== do overlaps found at m=800 survive at m=1600? ==")
    V8, V16 = load(800), load(1600)
    hits8 = run(V8, V8, 800, 800, thresh=0.15)
    hits8.sort()
    chosen = []
    for h in hits8:
        if all(abs(h[1] - g[1]) > 25 or abs(h[3] - g[3]) > 25 for g in chosen):
            chosen.append(h)
        if len(chosen) >= 8:
            break
    refine_rows = []
    for h in chosen:
        r8 = so.refine(V8, h)
        a, L, c, K = r8["a"], r8["L"], r8["c"], r8["K"]
        ratio = (len(V16) - 1) / (len(V8) - 1)
        a16, L16, c16, K16 = a * ratio, L * ratio, c * ratio, K * ratio
        t = np.linspace(0, 1, 201)
        def score(VBB, a_, L_, c_, K_):
            w = so.sample(VBB, a_ + t * L_)
            q = so.sample(VBB, c_ + t * K_)
            s, res, size = so.shape_score(w, q)
            return res
        r16 = score(V16, a16, L16, c16, K16)
        r16s = score(shuffled(V16, 1), a16, L16, c16, K16)
        # where does that land among random placements of the same lengths at m=1600?
        rng = np.random.default_rng(0)
        Lr = int(round(L16)); Kr = int(round(K16))
        vals = []
        for _ in range(300):
            a0 = rng.uniform(0, len(V16) - 1 - Lr)
            c0 = rng.uniform(0, len(V16) - 1 - Kr)
            vals.append(score(V16, a0, Lr, c0, Kr))
        pct = float(np.mean(np.array(vals) < r16))
        refine_rows.append(dict(resid8=r8["residual"], resid16=r16, resid16_shuffled=r16s,
                                percentile16=pct, a=a, L=L, c=c, K=K,
                                size=r8["size"], arg=r8["argC"]))
        print(f"  m=800 [{a:6.0f},{a+L:6.0f}] [{c:6.0f},{c+K:6.0f}] resid {r8['residual']:.4f}"
              f"  ->  m=1600 resid {r16:.4f} (shuffled {r16s:.3f}, "
              f"percentile {100*pct:.1f}%)")
    out["nesting"] = refine_rows

    # ---- 3. cross-m: block of one curve against arcs of another
    print("\n== cross-m search: block from one curve, arc from another ==")
    for (mB, mA) in ((800, 1600), (400, 800), (800, 400)):
        VB, VA = load(mB), load(mA)
        hits = run(VB, VA, mB, mA)
        r = np.array([h[0] for h in hits]) if hits else np.array([])
        # the corresponding interval, where the two samplings describe the same piece
        ratio = (mA - 1) / (mB - 1)
        best_overall = r.min() if r.size else None
        nontriv = []
        for h in hits:
            _, a, L, c, K = h
            if abs(c - a * ratio) > 0.05 * (mA - 1) or abs(K - L * ratio) > 0.10 * L * ratio:
                nontriv.append(h)
        rn = np.array([h[0] for h in nontriv]) if nontriv else np.array([])
        # sanity check: the two samplings of the same piece must nearly coincide
        t = np.linspace(0, 1, 121)
        a0 = 0.3 * (mB - 1)
        L0 = 0.2 * (mB - 1)
        w = so.sample(VB, a0 + t * L0)
        q = so.sample(VA, (a0 + t * L0) * ratio)
        _, resid_corr, _ = so.shape_score(w, q)
        sh = shuffled(VA, 1)
        hits_sh = coarse2(VB, sh, mB, mA)
        rsh = np.array([h[0] for h in hits_sh]) if hits_sh else np.array([])
        print(f"  {mB}->{mA}: best {best_overall if best_overall else float('nan'):.4f}"
              f"   best away from the correspondence {rn.min() if rn.size else float('nan'):.4f}"
              f"   (#<0.15 away: {int((rn < 0.15).sum()) if rn.size else 0})"
              f"   [corresponding piece scores {resid_corr:.4f}; shuffled arcs best "
              f"{rsh.min() if rsh.size else float('nan'):.3f}]")
        out[f"cross_{mB}_{mA}"] = dict(best=best_overall,
                                       best_nontrivial=float(rn.min()) if rn.size else None,
                                       n15_nontrivial=int((rn < 0.15).sum()) if rn.size else 0)

    Path(ROOT / "results/overlaps-vs-m.json").write_text(json.dumps(out, indent=1) + "\n")
    print("\nwrote results/overlaps-vs-m.json")


if __name__ == "__main__":
    main()
