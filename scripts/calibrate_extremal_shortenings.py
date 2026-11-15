"""Calibrate the shape test on an exact Koch quarter-shortening."""

from pathlib import Path

import numpy as np
import test_extremal_shortenings as check

ROOT = Path(__file__).resolve().parents[1]
ALPHA = np.log(3) / np.log(4)


def koch(x, depth=50):
    """Evaluate the usual four-map Koch arc by iterating its base-4 digits."""
    directions = np.array([
        1 + 0j, np.exp(1j * np.pi / 3),
        np.exp(-1j * np.pi / 3), 1 + 0j,
    ])
    factors = directions / 3
    offsets = np.r_[0, np.cumsum(directions)] / 3
    x = np.atleast_1d(np.asarray(x, dtype=float))
    y = x.copy()
    digits = np.zeros((depth, len(x)), dtype=int)
    for level in range(depth):
        y *= 4
        digit = np.minimum(y.astype(int), 3)
        digits[level] = digit
        y -= digit
    result = np.zeros(len(x), dtype=complex)
    for level in range(depth - 1, -1, -1):
        result = factors[digits[level]] * result + offsets[digits[level]]
    return result


def main():
    curve = koch(np.linspace(0, 1, 65))
    check.V = curve
    check.N = len(curve)
    check.DMATRIX = np.abs(curve[:, None] - curve[None, :])
    # A quarter at this grid is [0,16], since (k/16)/4 = k/64.
    _, best, _ = check.shape_error(0, 16)
    assert best["variant"] == "identity"
    assert best["sup"] < 1e-13, best

    quarter_ratio = abs(curve[16] - curve[0]) / (16 / 64) ** ALPHA
    lower, upper = 2 * 6 ** -ALPHA, 2 * np.sqrt(211) * 69 ** -ALPHA
    assert abs(quarter_ratio - 1) < 1e-14
    lines = [
        "Calibration of the shape comparison",
        "Exact quarter-copy of the Koch arc (65 nodes, recursive depth 50):",
        f"  identity match sup = {best['sup']:.17g}",
        f"  identity match RMS = {best['rms']:.17g}",
        "  PASS (< 1e-13)",
        "",
        f"The quarter-copy ratio r(0,1/4) is {quarter_ratio:.15f}. "
        "The established Koch global ratio range is",
        "[2*6^(-alpha), 2*sqrt(211)*69^(-alpha)] = "
        f"[{lower:.15f}, {upper:.15f}].",
        "Therefore a true self-shortening need not be a near-extremal chord.",
    ]
    text = "\n".join(lines) + "\n"
    print(text, end="")
    (ROOT / "results/extremal-shortenings-calibration.out").write_text(text)


if __name__ == "__main__":
    main()
