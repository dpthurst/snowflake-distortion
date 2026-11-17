"""Test whether long extremal chords select near-self-similar subarcs.

For a block v[a:b+1], the proposed shortening is compared against the full
saved configuration at matching parameter fractions (linear-in-index
interpolation only for this shape comparison). Both indices and conjugation
are reversed as four alternatives. All quantities use the original index grid.

Only two representatives of each extremal-chord class are selected. Nearby,
nearly duplicate endpoint pairs are suppressed. For each selected block we
also test ALL intervals with exactly the same index gap, rather than comparing
only with a random or arbitrary baseline.
"""

import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist, squareform

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/inputs/s_m800_plain.npy"
V = np.load(SOURCE)
assert V.shape == (800,) and np.isfinite(V).all()
V = (V - V[0]) / (V[-1] - V[0])
N = len(V)
ALPHA = np.log(3) / np.log(4)
II, JJ = np.triu_indices(N, 1)
RATIOS = np.abs(V[II] - V[JJ]) / (JJ - II) ** ALPHA
DMATRIX = squareform(pdist(np.column_stack([V.real, V.imag])))
VARIANTS = ("identity", "conjugation", "reversal", "reversal+conjugation")


def full_at(t):
    """Reference polyline in index parameter, only for comparing shapes."""
    q = np.asarray(t) * (N - 1)
    return np.interp(q, np.arange(N), V.real) + 1j * np.interp(
        q, np.arange(N), V.imag
    )


def transform(w, variant):
    if variant == 0:
        return w
    if variant == 1:
        return np.conj(w)
    if variant == 2:
        return 1 - w[::-1]
    return 1 - np.conj(w[::-1])


def shape_error(a, b):
    """Endpoint-gauged shape residual, normalized by the block's diameter."""
    assert 0 <= a < b < N
    w = (V[a:b + 1] - V[a]) / (V[b] - V[a])
    dia = DMATRIX[a:b + 1, a:b + 1].max() / abs(V[b] - V[a])
    reference = full_at(np.linspace(0, 1, len(w)))
    fits = []
    for variant in range(4):
        e = np.abs(w - transform(reference, variant))
        fits.append({
            "variant": VARIANTS[variant],
            "sup": float(e.max() / dia),
            "rms": float(np.sqrt(np.mean(e ** 2)) / dia),
        })
    best_sup = min(range(4), key=lambda k: (fits[k]["sup"], fits[k]["rms"]))
    best_rms = min(range(4), key=lambda k: (fits[k]["rms"], fits[k]["sup"]))
    return fits, fits[best_sup], fits[best_rms]


def representatives(cutoff, kind, count=2):
    """The top non-nearly-duplicate extremal chords in a separation band."""
    admissible = np.flatnonzero(JJ - II >= cutoff * (N - 1))
    if kind == "max":
        order = admissible[np.argsort(RATIOS[admissible])[::-1]]
        band_extreme = RATIOS[order[0]]
    else:
        order = admissible[np.argsort(RATIOS[admissible])]
        band_extreme = RATIOS[order[0]]
    pairs = []
    for k in order:
        a, b = int(II[k]), int(JJ[k])
        if all(max(abs(a - c), abs(b - d)) >= 12 for c, d in pairs):
            pairs.append((a, b))
        if len(pairs) == count:
            break
    return [{
        "a": a,
        "b": b,
        "ratio": float(RATIOS[np.flatnonzero((II == a) & (JJ == b))[0]]),
        "band_extreme_ratio": float(band_extreme),
        "relative_gap_to_band_extreme": float(abs(
            RATIOS[np.flatnonzero((II == a) & (JJ == b))[0]] - band_extreme
        ) / band_extreme),
    } for a, b in pairs]


def same_length_comparison(a, b):
    """Compare with every translation of [a,b] having the same index gap."""
    rows = [shape_error(c, c + b - a)[1:] for c in range(N - (b - a))]
    sup_rows = [row[0] for row in rows]
    rms_rows = [row[1] for row in rows]
    sup = np.array([r["sup"] for r in sup_rows])
    rms = np.array([r["rms"] for r in rms_rows])
    return {
        "number_of_intervals": len(rows),
        "minimum_sup": float(sup.min()),
        "best_sup_range_inclusive": [int(sup.argmin()), int(sup.argmin() + b - a)],
        "minimum_rms": float(rms.min()),
        "best_rms_range_inclusive": [int(rms.argmin()), int(rms.argmin() + b - a)],
        "median_sup": float(np.median(sup)),
        "sup_rank_lowest_is_one": int(np.sum(sup < sup_rows[a]["sup"]) + 1),
        "rms_rank_lowest_is_one": int(np.sum(rms < rms_rows[a]["rms"]) + 1),
    }


def main():
    output = {"source": str(SOURCE.relative_to(ROOT)), "alpha": float(ALPHA),
              "shape_comparison_interpolant": "linear in original index; shape test only",
              "near_extremal_band_tolerance": 0.02,
              "candidates": [], "near_extremal_bands": []}
    for cutoff in (0.2, 0.5):
        for kind in ("max", "min"):
            for rep in representatives(cutoff, kind):
                a, b = rep["a"], rep["b"]
                fits, best_sup, best_rms = shape_error(a, b)
                i, j = np.triu_indices(b - a + 1, 1)
                w = (V[a:b + 1] - V[a]) / (V[b] - V[a])
                r = np.abs(w[i] - w[j]) / (j - i) ** ALPHA
                row = {**rep, "cutoff": cutoff, "kind": kind,
                       "relative_length": float((b - a) / (N - 1)),
                       "similarity_scale": float(abs(V[b] - V[a])),
                       "similarity_angle_degrees": float(np.angle(V[b] - V[a], deg=True)),
                       "block_D": float(r.max() / r.min()),
                       "best_sup": best_sup, "best_rms": best_rms,
                       "all_fits": fits,
                       "same_length_baseline": same_length_comparison(a, b)}
                output["candidates"].append(row)
                print(json.dumps(row, indent=2), flush=True)
    # Widen the test: all chords within 2% of the extreme ratio in each band.
    # Report the best shape matches, not just the most extreme endpoint pairs.
    for cutoff in (0.2, 0.5):
        for kind in ("max", "min"):
            allowed = np.flatnonzero(JJ - II >= cutoff * (N - 1))
            ratios = RATIOS[allowed]
            selected = allowed[ratios >= ratios.max() / 1.02] if kind == "max" else \
                allowed[ratios <= ratios.min() * 1.02]
            scored = []
            for k in selected:
                a, b = int(II[k]), int(JJ[k])
                fits, best_sup, best_rms = shape_error(a, b)
                scored.append({"a": a, "b": b,
                               "ratio": float(RATIOS[k]), "fits": fits,
                               "best_sup": best_sup, "best_rms": best_rms})
            result = {"cutoff": cutoff, "kind": kind,
                      "number_of_chords_within_two_percent": len(scored)}
            for label, criterion in (("best_sup", "sup"), ("best_rms", "rms")):
                best = min(scored, key=lambda q: (q[label][criterion],
                                                 q[label]["rms" if criterion == "sup" else "sup"]))
                a, b = best["a"], best["b"]
                i, j = np.triu_indices(b - a + 1, 1)
                w = (V[a:b + 1] - V[a]) / (V[b] - V[a])
                r = np.abs(w[i] - w[j]) / (j - i) ** ALPHA
                best["block_D"] = float(r.max() / r.min())
                best["same_length_baseline"] = same_length_comparison(a, b)
                result[label] = best
            output["near_extremal_bands"].append(result)
            print(json.dumps(result, indent=2), flush=True)
    text = json.dumps(output, indent=2) + "\n"
    (ROOT / "results/extremal_shortenings.json").write_text(text)
    (ROOT / "results/extremal_shortenings.out").write_text(text)
    summary = []
    for band in output["near_extremal_bands"]:
        c = band["best_sup"]
        baseline = c["same_length_baseline"]
        summary.append({
            "cutoff": band["cutoff"], "kind": band["kind"],
            "near_extremal_chords": band["number_of_chords_within_two_percent"],
            "interval": [c["a"], c["b"]], "ratio": c["ratio"],
            "variant": c["best_sup"]["variant"],
            "sup": c["best_sup"]["sup"], "rms": c["best_sup"]["rms"],
            "same_length_count": baseline["number_of_intervals"],
            "same_length_rank": baseline["sup_rank_lowest_is_one"],
            "same_length_best": baseline["minimum_sup"],
            "block_D": c["block_D"],
        })
    (ROOT / "results/extremal-shortenings-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
