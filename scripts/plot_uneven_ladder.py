"""How the bounded-kappa ladder climbs: D against m, and the increment per doubling."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COL = {1: "#1f3f6f", 2: "#3a6ea5", 4: "#2c9c4a", 8: "#e08a2e", 16: "#c8342a", 32: "#8b2fa0"}

# polished uniform ladder from the note (SLP-certified configurations)
UNIFORM_POLISHED = [(40, 1.3523099782), (100, 1.4362283999), (200, 1.4908552600),
                    (400, 1.5417920000), (800, 1.5900261898), (1600, 1.6355308997)]


def collect():
    data = {}
    rows = json.loads((ROOT / "results/uneven-bounded.json").read_text())
    for r in rows:
        data.setdefault(r["par"], []).append((r["m"], r["D"]))
    for name, path in (("u", "results/uneven-uniform-large.json"),
                       ("k4", "results/uneven-kappa4-large.json")):
        p = ROOT / path
        if p.exists():
            for r in json.loads(p.read_text()):
                key = 1 if name == "u" else r["kappa"]
                data.setdefault(key, []).append((r["m"], r["D"]))
    p = ROOT / "results/uneven-kappa4-m400.json"
    if p.exists():
        for r in json.loads(p.read_text()):
            data.setdefault(r["kappa"], []).append((r["m"], r["D"]))
    p = ROOT / "results/uneven-kappa4-m800.json"
    if p.exists():
        for r in json.loads(p.read_text()):
            data.setdefault(r["kappa"], []).append((r["m"], r["D"]))
    for k in data:
        data[k] = sorted(set(data[k]))
    return data


def increments(points):
    """dD per doubling of m, plotted at the geometric mean of the pair"""
    out = []
    for (m1, d1), (m2, d2) in zip(points[:-1], points[1:]):
        out.append((np.sqrt(m1 * m2), (d2 - d1) / np.log2(m2 / m1)))
    return np.array(out)


def iid_points():
    """(m, min, max) over the i.i.d. draws: small m from the plain runs, large m from the homotopy"""
    out = {}
    for path in ("results/uneven-random.json",):
        p = ROOT / path
        if p.exists():
            for r in json.loads(p.read_text()):
                m = r["m"]
                out.setdefault(m, []).append(r["D"])
    p = ROOT / "results/uneven-random-large.json"
    if p.exists():
        for r in json.loads(p.read_text()):
            out.setdefault(r["m"], []).append(r["D"])
    return sorted((m, min(v), max(v)) for m, v in out.items())


def main():
    data = collect()
    fig, (ax, bz) = plt.subplots(1, 2, figsize=(10.6, 4.0),
                                 gridspec_kw=dict(width_ratios=[1.15, 1]))
    for k in sorted(data):
        ms = [p[0] for p in data[k]]
        Ds = [p[1] for p in data[k]]
        ax.semilogx(ms, Ds, "o-", ms=3.6, lw=1.1, color=COL[k], label=f"$\\kappa={k}$")
    ax.semilogx([p[0] for p in UNIFORM_POLISHED], [p[1] for p in UNIFORM_POLISHED],
                "s--", ms=4.2, lw=1.1, color="0.3", label="uniform, polished")
    iid = iid_points()
    if iid:
        ms = np.array([p[0] for p in iid], float)
        lo = np.array([p[1] for p in iid])
        hi = np.array([p[2] for p in iid])
        ax.plot(ms, lo, "D--", ms=3.2, lw=0.7, color="#000000",
                label="i.i.d. spacings (best of draws)", zorder=5)
        ax.errorbar(ms, lo, yerr=[np.zeros_like(lo), hi - lo], fmt="none",
                    ecolor="#000000", elinewidth=0.8, capsize=2.0, zorder=4)
    ax.set_xlabel("$m$", fontsize=11)
    ax.set_ylabel("best $D_m$", fontsize=11)
    ax.grid(alpha=0.25)
    ax.set_title("bounded gap ratio $\\kappa$: $D$ against $m$", fontsize=10)
    ax.legend(fontsize=7.6, loc="lower right", ncol=2, handlelength=1.8, labelspacing=0.35)

    inc4 = increments([p for p in data[4] if p[0] >= 24])
    inc1 = increments(UNIFORM_POLISHED)
    inc32 = increments([p for p in data[32]])
    bz.semilogx(inc4[:, 0], inc4[:, 1], "o-", ms=3.6, lw=1.1, color=COL[4], label="$\\kappa=4$")
    bz.semilogx(inc1[:, 0], inc1[:, 1], "s--", ms=4.2, lw=1.1, color="0.3",
                label="uniform, polished")
    bz.set_xlabel("$m$", fontsize=11)
    bz.set_ylabel("$\\Delta D$ per doubling", fontsize=11)
    bz.grid(alpha=0.25)
    bz.set_title("the climb barely decelerates", fontsize=10)
    bz.set_ylim(0.04, 0.095)
    bz.legend(fontsize=7.8, loc="upper right")
    # guide: the power-law fit to the kappa=4 ladder (Dinf = 1.768, s = 0.273)
    mm = np.logspace(np.log10(24), np.log10(800), 50)
    bz.semilogx(mm, 1.19 * 0.273 * np.log(2) * mm ** (-0.273), ":", lw=0.9, color="0.5",
                label="fit $cm^{-0.27}$")
    fig.tight_layout()
    out = ROOT / "figures/uneven-ladder.png"
    fig.savefig(out, dpi=170, facecolor="white", bbox_inches="tight",
                pil_kwargs={"optimize": True})
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size/1024:.0f} kB)")
    for k in sorted(data):
        print(f"  kappa={k}: " + ", ".join(f"{m}:{d:.5f}" for m, d in data[k]))


if __name__ == "__main__":
    main()
