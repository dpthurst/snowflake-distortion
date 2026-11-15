"""Directly check the point/mirror symmetry of one shortening of S_m800.

This deliberately computes only the whole curve and nodes 25,...,739, using
literal vectorised operations on the saved coordinates.  It does not search for
an optimal shortening or certify a limiting curve.
"""

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / ".ultramath/smwhdqmht/m1/s_m800_plain.npy"
ALPHA = np.log(3) / np.log(4)
A, B = 25, 739


def gauge(v):
    """Translate and dilate a nondegenerate arc so its endpoints are 0 and 1."""
    v = np.asarray(v, dtype=complex)
    assert v[-1] != v[0]
    return (v - v[0]) / (v[-1] - v[0])


def statistics(v, original_start):
    """Compute all quantities directly; index gaps refer to the original grid."""
    w = gauge(v)
    differences = np.abs(w[:, None] - w[None, :])
    diameter = float(differences.max())
    i, j = np.triu_indices(len(w), 1)
    ratios = np.abs(w[i] - w[j]) / (j - i) ** ALPHA

    point_error = np.abs(w + w[::-1] - 1)
    mirror_error = np.abs(w - (1 - np.conj(w[::-1])))
    rmax, rmin = float(ratios.max()), float(ratios.min())
    return {
        "node_range_inclusive": [
            int(original_start), int(original_start + len(v) - 1)
        ],
        "parameter_interval": [
            float(original_start / (len(ALL) - 1)),
            float((original_start + len(v) - 1) / (len(ALL) - 1)),
        ],
        "relative_parameter_length": float((len(v) - 1) / (len(ALL) - 1)),
        "number_of_nodes": len(v),
        "diameter": diameter,
        "D": rmax / rmin,
        "rmin": rmin,
        "rmin_pair_original_indices": [
            int(original_start + i[ratios.argmin()]),
            int(original_start + j[ratios.argmin()]),
        ],
        "rmax": rmax,
        "rmax_pair_original_indices": [
            int(original_start + i[ratios.argmax()]),
            int(original_start + j[ratios.argmax()]),
        ],
        "point_sup": float(point_error.max() / diameter),
        "point_sup_original_index": int(original_start + point_error.argmax()),
        "point_rms": float(np.sqrt(np.mean(point_error ** 2)) / diameter),
        "mirror_sup": float(mirror_error.max() / diameter),
        "mirror_rms": float(np.sqrt(np.mean(mirror_error ** 2)) / diameter),
    }


ALL = np.load(SOURCE)
assert ALL.shape == (800,) and np.isfinite(ALL).all(), "unexpected source data"
assert 0 < A < B < len(ALL), "the shortening must be proper and nondegenerate"

result = {
    "source": str(SOURCE.relative_to(ROOT)),
    "alpha": float(ALPHA),
    "definitions": {
        "D": "max |w_i-w_j|/(j-i)^alpha divided by the corresponding min",
        "point_sup": "max |w_i+w_(N-1-i)-1| / diam(w)",
        "point_rms": "sqrt(mean |w_i+w_(N-1-i)-1|^2) / diam(w)",
        "mirror_sup": "max |w_i-(1-conj(w_(N-1-i)))| / diam(w)",
        "mirror_rms": "sqrt(mean |w_i-(1-conj(w_(N-1-i)))|^2) / diam(w)",
    },
    "whole": statistics(ALL, 0),
    "shortened": statistics(ALL[A:B + 1], A),
}
print(json.dumps(result, indent=2))
(ROOT / "results/shortening_symmetry.json").write_text(
    json.dumps(result, indent=2) + "\n"
)
