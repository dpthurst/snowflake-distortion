"""Check the browser geometry in ``app/snowflake-align.html`` against numpy.

What this does: extracts the core-math block and the embedded curve data from
the built page, runs them under node, and compares every number with
``scripts/align_core.py``.  It also parses both inline scripts with
``node --check``, and corrupts one line of the core to confirm the comparison
can fail.

What this does *not* do: there is no browser here, so drawing, dragging and
layout are untested.  This is a logic test, not a GUI test.

    python3 scripts/test_align_app.py
"""

import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import align_core as core  # noqa: E402

NODE = Path.home() / "node/bin/node"
HTML = ROOT / "app/snowflake-align.html"
SCRATCH = ROOT / ".scratch"

CASES = [
    # whole curve, block source, a, b
    ("S-branch m=800", "S-branch m=800", 455, 689),
    ("S-branch m=800", "S-branch m=800", 259, 511),
    ("S-branch m=800", "S-branch m=800", 165, 623),
    ("S-branch m=800", "S-branch m=800", 25, 739),
    ("S-branch m=800", "S-branch m=800", 0, 799),
    ("S-branch m=800", "S-branch m=800", 3, 4),
    ("S-branch m=100", "S-branch m=100", 0, 50),
    ("S-branch m=40", "S-branch m=40", 7, 33),
    ("S-branch m=400", "S-branch m=800", 186, 586),
    ("three-piece m=800", "three-piece m=800", 100, 700),
]


def extract(html):
    core_js = re.search(r"/\* ==== CORE MATH BEGIN ==== \*/(.*?)/\* ==== CORE MATH END ==== \*/",
                        html, re.S)
    curves = re.search(r"const CURVES = (\{.*?\});\n", html, re.S)
    if not core_js or not curves:
        raise SystemExit("could not extract the core or the data from the page")
    return core_js.group(1), json.loads(curves.group(1))


def script_bodies(html):
    return re.findall(r"<script>(.*?)</script>", html, re.S)


def run_node(core_js, curves, cases):
    SCRATCH.mkdir(exist_ok=True)
    exports = ("\nmodule.exports = {ALPHA, VARIANTS, gauge, interpAt, blockAt, wholeAt, knots,\n"
               "  supDiff, rmsDiff, residuals, procrustes, distortion, ladder, diamond,\n"
               "  stripParams, fitCopy, copyResidual, profileAt, fitTranslation,\n"
               "  bboxOf, setDistance, pointSegDist};\n")
    (SCRATCH / "align_core.js").write_text(core_js + exports)
    driver = """
const core = require('./align_core.js');
const fs = require('fs');
const { CURVES, CASES } = JSON.parse(fs.readFileSync('./align_cases.json', 'utf8'));
function arr(a) { return { x: Float64Array.from(a.x), y: Float64Array.from(a.y) }; }
const COPY_CASES = JSON.parse(fs.readFileSync('./align_cases.json', 'utf8')).COPY_CASES || [];
const outCopy = [];
for (const [w, s, a, b] of COPY_CASES) {
  const W = arr(CURVES[w]), S = arr(CURVES[s]);
  const fit = core.fitCopy(W.x, W.y, S.x, S.y, a, b);
  const Wg = core.gauge(W.x, W.y);
  const t = new Float64Array(1201); for (let i = 0; i < 1201; i++) t[i] = i / 1200;
  const wS = core.blockAt(S.x, S.y, a, b, t);
  const inside = core.copyResidual(Wg.x, Wg.y, { t, x: wS.x, y: wS.y },
                                   fit.cx, fit.cy, fit.dx, fit.dy, fit.s0, fit.s1);
  outCopy.push({ s0: fit.s0, s1: fit.s1, cx: fit.cx, cy: fit.cy, dx: fit.dx, dy: fit.dy,
                 sup: inside.sup, rms: inside.rms, natural_sup: fit.natural.sup });
}
const out = [];
for (const [w, s, a, b] of CASES) {
  const W = arr(CURVES[w]), S = arr(CURVES[s]);
  const res = core.residuals(W.x, W.y, S.x, S.y, a, b);
  let best = core.VARIANTS[0];
  for (const k of core.VARIANTS) if (res[k].sup < res[best].sup) best = k;
  const fit = core.procrustes(W.x, W.y, S.x, S.y, a, b, best);
  const seg = { x: S.x.slice(a, b + 1), y: S.y.slice(a, b + 1) };
  const rho = Math.hypot(S.x[b] - S.x[a], S.y[b] - S.y[a]);
  const lad = core.ladder(S.x, S.y, b - a);
  const r_ab = rho / Math.pow(b - a, core.ALPHA);
  let below = 0; for (let i = 0; i < lad.length; i++) if (i !== a && lad[i] < r_ab) below++;
  const pct = lad.length > 1 ? below / (lad.length - 1) : null;
  /* the residual strip must be drawn on a set containing both knot sets, so that
     its plotted maximum is the exact sup and not a sampled one */
  const tp = core.stripParams(W.x.length, b - a);
  const Wg = core.gauge(W.x, W.y);
  const gs = core.wholeAt(Wg.x, Wg.y, tp, best);
  const ws = core.blockAt(S.x, S.y, a, b, tp);
  let smax = 0;
  for (let i = 0; i < tp.length; i++) smax = Math.max(smax, Math.hypot(ws.x[i] - gs.x[i], ws.y[i] - gs.y[i]));
  /* symmetric distance between the two polylines, as the page computes it */
  const NS = 601, tsd = new Float64Array(NS);
  for (let i = 0; i < NS; i++) tsd[i] = i / (NS - 1);
  const wsd = core.blockAt(S.x, S.y, a, b, tsd);
  const Wgd = core.gauge(W.x, W.y);
  const gsd = core.wholeAt(Wgd.x, Wgd.y, tsd, best);
  const tv = new Float64Array(W.x.length);
  for (let k = 0; k < tv.length; k++) tv[k] = k / (W.x.length - 1);
  const gv = core.wholeAt(Wgd.x, Wgd.y, tv, best);
  const bt = new Float64Array(b - a + 1);
  for (let k = 0; k <= b - a; k++) bt[k] = k / (b - a);
  const bv = core.blockAt(S.x, S.y, a, b, bt);
  const boxG = core.bboxOf(gv.x, gv.y, 0, gv.x.length - 1);
  const boxB = core.bboxOf(bv.x, bv.y, 0, bv.x.length - 1);
  const sdF = core.setDistance(wsd.x, wsd.y, gv.x, gv.y, boxG).sup;
  const sdB = core.setDistance(gsd.x, gsd.y, bv.x, bv.y, boxB).sup;
  out.push({ setdist: [sdF, sdB, Math.max(sdF, sdB)],
             strip_max: smax,
             sup: Object.fromEntries(core.VARIANTS.map(k => [k, res[k].sup])),
             rms: Object.fromEntries(core.VARIANTS.map(k => [k, res[k].rms])),
             best, free_sup: fit.sup, absC: fit.absC, argC: fit.argC,
             D_block: core.distortion(seg.x, seg.y),
             s_block: core.diamond(S.x, S.y, a, b) / rho,
             pct, r_ab });
}
process.stdout.write(JSON.stringify({ cases: out, copy: outCopy }));
"""
    (SCRATCH / "align_run.js").write_text(driver)
    cases_in = {"CURVES": curves, "CASES": cases, "COPY_CASES": COPY_CASES}
    (SCRATCH / "align_cases.json").write_text(json.dumps(cases_in))
    proc = subprocess.run([str(NODE), str(SCRATCH / "align_run.js")],
                          cwd=SCRATCH, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit("node failed:\n" + proc.stderr[-4000:])
    return json.loads(proc.stdout)


def python_side(curves, case):
    w, s, a, b = case
    W = np.array(curves[w]["x"]) + 1j * np.array(curves[w]["y"])
    S = np.array(curves[s]["x"]) + 1j * np.array(curves[s]["y"])
    src = None if s == w else S
    res = core.residuals(W, a, b, src)
    best = min(core.VARIANTS, key=lambda k: res[k]["sup"])
    c, free_sup = core.procrustes(W, a, b, best, src)
    seg = S[a:b + 1]
    rho = abs(S[b] - S[a])
    lad = core.ladder(S, b - a)
    r_ab = rho / (b - a) ** core.ALPHA
    sd = core.set_distance_sup(W, a, b, src, kind=best)
    return dict(setdist=list(sd), strip_max=res[best]["sup"], sup={k: res[k]["sup"] for k in core.VARIANTS},
                rms={k: res[k]["rms"] for k in core.VARIANTS},
                best=best, free_sup=free_sup, absC=abs(c), argC=float(np.angle(c)),
                D_block=core.distortion(seg),
                s_block=float(np.abs(seg[:, None] - seg[None, :]).max() / rho),
                pct=(float((lad[np.arange(lad.size) != a] < r_ab).mean())
                     if lad.size > 1 else None), r_ab=r_ab)


def compare(js, py, case, tol=1e-9, report=None):
    bad = []
    for key in ("sup", "rms"):
        for k in js[key]:
            a, b = js[key][k], py[key][k]
            if abs(a - b) > tol * max(1.0, abs(b)):
                bad.append(f"{key}[{k}] js={a!r} py={b!r}")
    for i, v in enumerate(js.get("setdist", [])):
        if abs(v - py["setdist"][i]) > tol * max(1.0, abs(py["setdist"][i])):
            bad.append(f"setdist[{i}] js={v!r} py={py['setdist'][i]!r}")
    for key in ("strip_max", "free_sup", "absC", "D_block", "s_block", "pct", "r_ab"):
        if js[key] is None and py[key] is None:
            continue
        if js[key] is None or py[key] is None:
            bad.append(f"{key} js={js[key]!r} py={py[key]!r} (one side undefined)")
            continue
        if abs(js[key] - py[key]) > tol * max(1.0, abs(py[key])):
            bad.append(f"{key} js={js[key]!r} py={py[key]!r}")
    if js["best"] != py["best"]:
        bad.append(f"best js={js['best']} py={py['best']}")
    if abs(js["argC"] - py["argC"]) > 1e-9 and abs(abs(js["argC"] - py["argC"]) - 2 * np.pi) > 1e-9:
        bad.append(f"argC js={js['argC']} py={py['argC']}")
    if bad:
        msg = f"case {case}: " + "; ".join(bad)
        (report.append(msg) if report is not None else (_ for _ in ()).throw(SystemExit(msg)))
    return not bad


def build_copy_cases(curves):
    """Synthetic curves that are exact (or only chord-approximate) similarity
    images of a window of the m=800 curve, under a known affine parameter."""
    import copy as _copy

    W = np.array(curves["S-branch m=800"]["x"]) + 1j * np.array(curves["S-branch m=800"]["y"])
    Wg = core.gauge(W)
    m = len(W)
    c, d = 0.37 * np.exp(-1.2j), 0.21 + 0.33j
    out = _copy.deepcopy(curves)
    cases = []
    for name, (j0, q, N) in {"synthetic exact": (100, 1, 200),
                             "synthetic coarse": (100, 2, 200)}.items():
        nodes = c * Wg[j0 + np.arange(N) * q] + d
        out[name] = dict(x=[round(float(z.real), 12) for z in nodes],
                         y=[round(float(z.imag), 12) for z in nodes],
                         D=core.distortion(nodes), src="synthetic")
        cases.append(("S-branch m=800", name, 0, N - 1, j0 / (m - 1),
                      q * (N - 1) / (m - 1)))
    cases.append(("S-branch m=800", "S-branch m=40", 0, 39, None, None))
    return out, cases


def main():
    html = HTML.read_text()
    core_js, curves = extract(html)
    curves, raw_copy_cases = build_copy_cases(curves)
    global COPY_CASES, COPY_TRUTH
    COPY_CASES = [c[:4] for c in raw_copy_cases]
    COPY_TRUTH = {c[1]: (c[4], c[5]) for c in raw_copy_cases if c[4] is not None}

    for i, body in enumerate(script_bodies(html)):
        tmp = SCRATCH / f"app_script_{i}.js"
        SCRATCH.mkdir(exist_ok=True)
        tmp.write_text(body)
        chk = subprocess.run([str(NODE), "--check", str(tmp)], capture_output=True, text=True)
        if chk.returncode != 0:
            raise SystemExit(f"script block {i} does not parse:\n{chk.stderr[-3000:]}")
    print(f"both inline scripts parse ({len(script_bodies(html))} blocks)")

    js = run_node(core_js, curves, CASES)["cases"]
    failures = []
    for case, j in zip(CASES, js):
        p = python_side(curves, case)
        ok = compare(j, p, case, report=failures)
        print(f"  {'ok ' if ok else 'BAD'} whole={case[0]:18s} source={case[1]:18s} "
              f"[{case[2]},{case[3]}]  sup={p['sup'][p['best']]:.9f} best={p['best']}")
    if failures:
        raise SystemExit("mismatches:\n  " + "\n  ".join(failures))

    # ---- copy mode: similarity + affine reparametrisation
    js_copy = run_node(core_js, curves, CASES)["copy"]
    for case, j in zip(COPY_CASES, js_copy):
        w, src, a, b = case
        if src in COPY_TRUTH:
            j["s0_true"], j["s1_true"] = COPY_TRUTH[src]
        W = np.array(curves[w]["x"]) + 1j * np.array(curves[w]["y"])
        S = np.array(curves[src]["x"]) + 1j * np.array(curves[src]["y"])
        t = np.linspace(0, 1, 1201)
        sup, rms = core.copy_residual(W, a, b, S, j["s0"], j["s1"],
                                      complex(j["cx"], j["cy"]), complex(j["dx"], j["dy"]), t=t)
        ref = core.fit_copy_multires(W, a, b, S)
        tag = f"{w} vs {src} [{a},{b}]"
        assert abs(sup - j["sup"]) < 1e-9, f"{tag}: python recompute {sup} vs js {j['sup']}"
        if src.startswith("synthetic exact"):
            assert j["sup"] < 1e-6, f"{tag}: exact copy not recovered, sup={j['sup']}"
            assert abs(j["s0"] - j["s0_true"]) < 2e-3 and abs(j["s1"] - j["s1_true"]) < 2e-3, \
                f"{tag}: parameters s0={j['s0']} s1={j['s1']} not recovered"
            print(f"  ok  copy fit recovers an exact copy: {tag}  sup={j['sup']:.2e} "
                  f"s=({j['s0']:.4f},{j['s1']:.4f}) truth ({j['s0_true']:.4f},{j['s1_true']:.4f})")
        else:
            assert j["sup"] <= j["natural_sup"] + 1e-9, \
                f"{tag}: fit {j['sup']} is worse than the index correspondence {j['natural_sup']}"
            assert j["sup"] <= ref["sup"] + 0.02, \
                f"{tag}: js {j['sup']} much worse than the reference {ref['sup']}"
            print(f"  ok  copy fit: {tag}  sup={j['sup']:.5f} (index correspondence "
                  f"{j['natural_sup']:.5f}, reference {ref['sup']:.5f})")

    # the comparison must be able to fail: flip a sign in the free-fit multiplier
    broken = core_js.replace("numY += ax * by - ay * bx;", "numY += ay * bx - ax * by;")
    assert broken != core_js, "corruption target not found"
    js_bad = run_node(broken, curves, CASES[:3])["cases"]
    caught = []
    for case, j in zip(CASES[:3], js_bad):
        compare(j, python_side(curves, case), case, report=caught)
    print(f"corruption test: {len(caught)} mismatches detected after flipping the sign in the "
          f"free-fit numerator -> {'PASS' if caught else 'FAIL'}")
    if not caught:
        raise SystemExit("the comparison did not notice a deliberate sign flip")

    # and the strip must not silently fall back to a uniform grid
    coarse = core_js.replace("const all = Array.from(k).concat(Array.from(u)).sort((p, q) => p - q);",
                             "const all = Array.from(u).sort((p, q) => p - q);")
    assert coarse != core_js, "strip corruption target not found"
    js_coarse = run_node(coarse, curves, CASES[:4])["cases"]
    caught2 = []
    for case, j in zip(CASES[:4], js_coarse):
        compare(j, python_side(curves, case), case, report=caught2)
    print(f"strip test: {len(caught2)} mismatches when the knots are dropped from the strip -> "
          f"{'PASS' if caught2 else 'FAIL'}")
    if not caught2:
        raise SystemExit("a uniform-only strip went unnoticed")
    # the distance metric must notice a truncated polyline search
    one_seg = core_js.replace("for (let i = 0; i < boxes.length; i++) {\n    const b = boxes[i];",
                              "for (let i = 0; i < 1; i++) {\n    const b = boxes[i];")
    if one_seg == core_js:
        one_seg = core_js.replace("for (let i = 0; i < boxes.length; i++) {", "for (let i = 0; i < 1; i++) {")
    assert one_seg != core_js, "set-distance corruption target not found"
    js_one = run_node(one_seg, curves, CASES)["cases"]
    caught_sd = []
    for case, j in zip(CASES, js_one):
        compare(j, python_side(curves, case), case, report=caught_sd)
    print(f"set-distance test: {len(caught_sd)} mismatches with only the first segment "
          f"searched -> {'PASS' if caught_sd else 'FAIL'}")
    if not caught_sd:
        raise SystemExit("a truncated polyline search went unnoticed")

    # and the copy fit must notice a broken parameter correspondence
    bad_sigma = core_js.replace("specific", "specific")
    bad_sigma = core_js.replace("for (let i = 0; i < n; i++) ss[i] = s0 + s1 * w.t[i];",
                                "for (let i = 0; i < n; i++) ss[i] = s0 + w.t[i];")
    assert bad_sigma != core_js, "sigma corruption target not found"
    js_bad_sigma = run_node(bad_sigma, curves, CASES)["copy"]
    broke = False
    for case, j in zip(COPY_CASES, js_bad_sigma):
        if case[1] == "synthetic exact" and j["sup"] > 1e-6:
            broke = True
    print(f"copy test: exact recovery lost when the affine parameter is broken -> "
          f"{'PASS' if broke else 'FAIL'}")
    if not broke:
        raise SystemExit("a broken affine correspondence was not noticed")
    print("all comparisons agree")


if __name__ == "__main__":
    main()
