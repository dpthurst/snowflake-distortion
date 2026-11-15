"""Build ``app/snowflake-align.html``: a self-contained alignment explorer.

The page embeds the saved configurations (small ``m`` ladder plus the arc and
three-piece branches) and all geometry is computed in the browser, so the file
opens offline with no server.  Rebuild after changing the template or the data:

    python3 scripts/build_align_app.py
"""

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys_path = ROOT / "scripts"

import sys
sys.path.insert(0, str(sys_path))
import align_core as core  # noqa: E402

U = ROOT / ".ultramath/smwhdqmht/m1"
M8 = ROOT / ".ultramath/smwh6zmc9/m8"
SCRATCH = ROOT / ".scratch"

SOURCES = {
    "S-branch m=40": (SCRATCH / "final_m40.npz", "W", 0),
    "S-branch m=100": (SCRATCH / "final_m100.npz", "W", 0),
    "S-branch m=150": (SCRATCH / "cert150.npy", None, None),
    "S-branch m=200": (SCRATCH / "descended_m200.npy", None, None),
    "S-branch m=400": (U / "s_m400.npy", None, None),
    "S-branch m=800": (U / "s_m800_plain.npy", None, None),
    "S-branch m=1600": (U / "s_m1600_plain.npy", None, None),
    "arc branch m=400": (U / "arc_m400.npy", None, None),
    "arc branch m=800": (U / "arc_m800.npy", None, None),
    "arc branch m=1600": (U / "arc_m1600.npy", None, None),
    "three-piece m=400": (M8 / "cand_m400.npy", None, None),
    "three-piece m=800": (M8 / "cand_m800.npy", None, None),
    "three-piece m=1600": (M8 / "cand_m1600.npy", None, None),
    "smoothed S m=800": (U / "S_m800_smooth.npy", None, None),
}

PRESETS = {
    "S-branch m=800": [
        [259, 511, "best residual, rel. length 0.32"],
        [165, 623, "best polyline residual below 0.75"],
        [277, 479, "near-maximal chord, best of its length"],
        [112, 669, "near-maximal chord, long"],
        [455, 689, "maximal chord (L>=160)"],
        [356, 607, "minimal chord (L>=160)"],
        [190, 590, "best RMS at length 400"],
        [25, 739, "the earlier symmetry check"],
    ],
    "S-branch m=1600": [
        [347, 1146, "best cross-resolution window (m=1600 to m=800)"],
        [358, 1157, "maximal chord, half length"],
        [666, 1465, "minimal chord, half length"],
    ],
    "S-branch m=400": [
        [0, 200, "first half"],
        [0, 100, "first quarter"],
        [100, 300, "middle half"],
    ],
    "S-branch m=100": [
        [0, 50, "first half"],
        [0, 25, "first quarter"],
    ],
}


def load(path, field, index):
    if path.suffix == ".npz":
        z = np.load(path)
        return np.asarray(z[field][index], dtype=complex)
    return np.asarray(np.load(path), dtype=complex)


def main():
    curves = {}
    for name, (path, field, index) in SOURCES.items():
        if not path.exists():
            print(f"skip {name}: missing {path.name}")
            continue
        v = core.gauge(load(path, field, index))
        curves[name] = {
            "x": [round(float(a), 7) for a in v.real],
            "y": [round(float(b), 7) for b in v.imag],
            "D": core.distortion(v),
            "src": str(path.relative_to(ROOT)),
        }
    data = json.dumps({"x": 0}, separators=(",", ":"))  # placeholder for lint

    packed = {
        name: {"x": c["x"], "y": c["y"], "D": c["D"], "src": c["src"]}
        for name, c in curves.items()
    }
    template = (ROOT / "scripts/align_app_template.html").read_text()
    html = (template
            .replace("__CURVES__", json.dumps(packed, separators=(",", ":")))
            .replace("__PRESETS__", json.dumps(PRESETS, separators=(",", ":"))))
    out = ROOT / "app/snowflake-align.html"
    out.parent.mkdir(exist_ok=True)
    out.write_text(html)
    print(f"wrote {out.relative_to(ROOT)}  ({out.stat().st_size/1024:.0f} kB, {len(curves)} curves)")
    for name, c in curves.items():
        print(f"  {name:22s} m={len(c['x']):5d}  D={c['D']:.9f}")


if __name__ == "__main__":
    main()
