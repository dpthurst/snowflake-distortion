# The α-snowflake of an interval: distortion, the Koch extremal, and the finite-$m$ ladder

For $\alpha=\log 3/\log 4$, the bi-Hölder distortion of an arc $\gamma\colon[0,1]\to\mathbb R^2$ is
$$D(\gamma)=\frac{\sup_{x\ne y}|\gamma(x)-\gamma(y)|/|x-y|^\alpha}{\inf_{x\ne y}|\gamma(x)-\gamma(y)|/|x-y|^\alpha},$$
and $D_\infty$ is its infimum over injective arcs. This repository contains two notes and the code
and data behind them.

**Established.** The Koch parametrisation has $D=\sqrt{633}\,23^{-\alpha}=2.096805160875$, with both
constants sharp (`notes/koch-extremals.pdf`, computer-assisted, with the endpoint value formalised in
Lean). Two reversible self-similar arcs sit below it at $2.042$ and $2.0179$ (measured, not
certified). The first rungs of the finite-$m$ ladder are certified globally by interval
branch-and-bound: $D_4=1.062765923076\ldots$, $D_5\in[1.093090076,1.093091076]$,
$D_6\in[1.127430051,1.128130050]$, hence $1.1274\le D_\infty\le2.0968$. Best known rungs climb from
$1.352$ at $m=40$ to $1.636$ at $m=1600$.

**Open.** Any better lower bound (a second-order bound might reach $m=7,8$); whether either explicit
arc is optimal; which macro-shape, if any, the ladder converges to; and the critical exponent
$\alpha=\tfrac12$, where a finite-distortion embedding of $[0,1]$ with distance $|x-y|^{1/2}$ is
neither known nor excluded.

## Layout

| path | contents |
|---|---|
| `notes/improving-koch-extremals.{tex,pdf}` | the main note: explicit arcs, the ladder, the certified rungs, heat maps, repeated blocks |
| `notes/koch-extremals.{tex,pdf}` | the companion: sharp constants for the Koch parametrisation |
| `certificate/` | the enclosures behind the companion's two constants, as cited by its §Certification: the lower-bound branch-and-bound, the maximum computation with its output, and the independent minimum-local checks |
| `figures/` | every figure of both notes, plus `figures/README.md` saying what each one is |
| `data/` | the node positions behind the figures and Table 1 (S branch, one-sided branch, $m=4,5,6$ optima), and `data/inputs/` for the other configurations the analysis scripts read |
| `results/` | working reports and JSON/`.npz` outputs |
| `scripts/` | the code: curve constructions, the minimax SLP, the certificates, the searches, the plotting |
| `lean/` | the Lean formalisation of the Koch endpoint value |
| `app/` | a browser tool for aligning blocks of a curve by hand (generated; see `scripts/build_align_app.py`) |

## Reproducing

The scripts are plain Python (numpy, scipy, mpmath) and run from the repository root, e.g.

```
python3 scripts/explicit_curves.py          # the two reversible arcs and the Koch arc
python3 scripts/plot_explicit_and_ladder.py # Figures 1, 7, 8, 9 from data/
python3 scripts/certify_mm.py --m 5 --eps 1e-6
python3 scripts/plot_heatmaps_selfsimilar.py
```

Lean: `lean/KochSnowflake.lean` compiles against Mathlib 4.35 (`lean/README.md`).
