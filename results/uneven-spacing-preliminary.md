# Uneven parameter spacings — preliminary

**Status.** Preliminary, deliberately not in either note (`notes/improving-koch-extremals.tex` had a
paragraph and figure on this and they were removed on 2026-11-16; the history still has them). To be
written up only once we agree what the facts are. Scripts: `scripts/optimise_uneven.py`,
`scripts/uneven_spacing.py`, `scripts/analyse_random_spacing.py`, `scripts/plot_uneven.py`,
`scripts/plot_uneven_random.py`. Data: `results/uneven-*.json`, `results/uneven_*.npz`.
Figures: `figures/uneven-spacing.png`, `figures/uneven-random-shapes.png`,
`figures/uneven-random-diagnostics.png`.

The problem with free parameters: $0=x_0<\dots<x_{m-1}=1$ and $v_i\in\mathbb{C}$, minimise the
spread of $r_{ij}=|v_i-v_j|/(x_j-x_i)^\alpha$ over all pairs.

## 1. Free spacings are degenerate

Three points already give $D=1$ exactly (a triangle with sides $c,c,2^\alpha c$ exists, with
$c=3^{-\alpha}$). Splitting them into near-duplicates — $\varepsilon$ apart in parameter and
$c\,\varepsilon^\alpha$ apart in space — keeps every ratio at 1, so the infimum is $1$ for every
$m\ge3$, approached only as the duplicates coalesce, and no curve is approximated. The optimiser
finds exactly this: $D=1.0008$ at $m=6$, gaps from $2.5\times10^{-5}$ to $0.56$, the six points in
three pairs. **The free problem is empty**; every statement below needs a constraint on the gaps.

## 2. Bounded gap ratio $\rho=\max\Delta x/\min\Delta x$

$\rho=1$ is the uniform ladder. For each fixed $\rho$ the compactness argument still gives
$D_m(\rho)\to D_\infty$, since comparable gaps keep the parameters asymptotically uniform.

```
          m=8      m=12     m=16     m=24     m=32     m=48     m=64     m=96
 rho=1   1.15661  1.22013  1.24535  1.29708  1.33192  1.37189  1.40031  1.43403
 rho=2   1.15012  1.20140  1.23488  1.28427  1.31850
 rho=4   1.11664  1.18004  1.21941  1.27031  1.30524
 rho=8   1.09427  1.15733  1.18935  1.25117  1.28663
 rho=16  1.08978  1.13516  1.17524  1.23142  1.27044
 rho=32  1.06277  1.09517  1.16348  1.20199  1.26423
```

Unevenness lowers every rung by a few per cent. The optimal spacings are *hierarchical*, not smooth:
at $\rho=32$, $m=32$ the gaps alternate between $0.0026$ and $0.082$ — clusters of near-duplicates
separated by jumps of the maximal size (`figures/uneven-spacing.png`, right panel). Imposing a smooth
profile instead is catastrophic: $x_i=(i/(m-1))^{1/2}$ gives $D=5.14$ at $m=12$ and
$x_i=(i/(m-1))^{2}$ gives $11.16$, against $1.22$ for uniform. It is the *structure* of the
unevenness that matters, not its size.

## 3. Random fixed spacings

$x_i$ = sorted uniform random, then frozen; only the positions are optimised. Three draws,
$m=16\dots64$.

```
   m     uniform    draw 0    draw 1    draw 2      gap ratio (draws)
  16     1.24535   1.21758   1.20414   1.21832      47 - 288
  24     1.29708   1.26336   1.26667   1.27676      48 - 216
  32     1.33192   1.29296   1.27874   1.30292     212 - 2199
  48     1.37189   1.34271   1.34781   1.34688     155 - 1934
  64     1.40031   1.36897   1.35468   1.37138      86 - 1319
```

Two things. **Random spacings beat the uniform grid at every $m$**, by 2–3%: the uniform
parametrisation is a poor sampling of the extremal curve, not the best $m$-point sampling. And the
values still increase with $m$, as they must ($D_m\le D_{m+1}\le D_\infty$).

**The shapes do not separate.** `figures/uneven-random-shapes.png` shows draw 0 at $m=16,32,64$
against the $m=800$ uniform configuration (grey). The coarse configuration is a polygon whose
vertices sit on the same S profile; by $m=64$ the points hug it. No second shape appears.

Quantitatively, with the conjugation ambiguity resolved by taking the better of $v$ and $\bar v$:

* max distance from a configuration to the $m=800$ polyline: random draws $0.045$–$0.11$
  (mean $0.013$–$0.053$); uniform rungs $0.034$–$0.12$ (mean $0.010$–$0.051$). The same range,
  falling with $m$ in both (`figures/uneven-random-diagnostics.png`, right).
* shape distances (max point-to-polyline, both directions) between the three $m=64$ draws:
  $0.068$, $0.069$, $0.081$. Between a draw and the uniform $m=64$: $0.054$, $0.081$, $0.085$.
  Between uniform rungs themselves ($m=96$ vs $m=800$): $0.054$.

So at these $m$ all the configurations, uniform and uneven, sit in one shape class to within
$0.05$–$0.08$ of the diameter: the same S profile, unevenly sampled and further along.

**A caution about the test.** "Qualitatively different from the exactly self-similar curves we have
looked at" would not by itself indicate a new type of optimum: the *uniform* S family is already not
self-similar in the measured sense (best sub-arc similarity residual $0.13$,
`results/shortening-similarity.md`). Differing from Koch, the two-piece and the three-piece arcs is a
property of the S branch already.

## 4. Caveats, and what would settle it

* All values are local optima from a few starts. At fixed $\rho$ they are not monotone in $m$
  ($\rho=32$: $1.0628$ at $m=8$, $1.0952$ at $m=12$), which cannot happen for exact optima. Seeding
  $m=12,16$ from the $m=8$ optimum improves them slightly ($1.1352\to1.1223$, $1.1635\to1.1516$) but
  not to $1.0628$.
* The shape test is limited by not knowing the limit curve: the uneven configurations are *better*
  than the $m=800$ uniform one, so they cannot coincide with it, and the distance to it is a crude
  proxy.
* What would decide it: run the random-spacing family to $m\approx200$–$400$ and see whether the
  pairwise shape distances keep falling ($0.07$ at $m=64$, against $0.054$ between uniform rungs at
  $m\ge96$); or compare sampling *densities* (the parameter measure) rather than point sets; or run
  the sub-arc self-similarity diagnostic directly on the uneven optima, which is the test that
  distinguishes "a different kind of optimum" from "a different sampling of the same curve".

## 5. Is the gain "just the prevalence of large arcs"?

Suggested hypothesis: at fixed $m$ an uneven parametrisation might win simply because it offers
long parameter spans whose ratios are slack, i.e. because "large arcs" are prevalent. The data says
no, on three counts (all at $m=32$ unless stated).

**The optimiser does not make large arcs.** At $\rho=32$ the longest gap is $2.5\times$ the mean gap
($0.082$ against $0.032$) and the shortest is $0.08\times$: it *clusters* points, it does not open
long arcs. Across $m=8\dots32$ at $\rho=32$ the largest gap is $2.4$–$4.8\times$ the mean.

**The gain does not need them.** Gaps confined to $[0.65,1.4]\times$ the mean ($\rho=2$, so no long
arcs at all) already recover $7$–$20\%$ of the $\rho=1\to32$ gain, $\rho=4$ recovers $28$–$43\%$ and
$\rho=16$ recovers $68$–$91\%$; the gain grows roughly linearly in $\log\rho$ with no threshold.

**The binding pairs are not the long ones.** At $m=32$:

```
                     D        max r at span   min r at span    r over spans<0.05   r over spans>0.5
 uniform          1.33192    0.387 (1.1x)     0.452 (1.3x)       0.789-0.938        0.789-1.050
 rho=2            1.31850    0.641 (1.8x)     0.535 (1.5x)       0.787-0.924        0.787-1.038
 rho=32           1.26423    0.805 (1.9x)     0.538 (1.2x)       0.825-1.043        0.825-1.043
 random draw 0    1.29296     -                -                  0.804-1.038        0.804-1.039
```

(spans in ordinary units, with "x" the span divided by the mean span of all pairs.)

The binding ratios sit at spans of $1.1$–$1.9\times$ the mean in every case --- never at the longest
arcs --- and the long-span ratios are if anything *narrower* than the rest, not slack in a way that
helps.

**What the gain actually is.** It is a rise of the *minimum* ratio: the smallest ratio goes from
$0.789$ (uniform) to $0.825$ ($\rho=32$), while the largest barely moves ($1.050\to1.043$). In the
uniform optimum the short-span pairs sit at the bottom of the ratio range ($0.789$--$0.938$ for spans
below $0.05$, against $1.050$ for long spans): the small-scale chords are systematically too short.
Uneven spacing lets the finite configuration place its points so that short and long spans have the
*same* ratio range ($0.825$--$1.043$ in both), i.e. it equalises across scales at fixed $m$, which the
uniform grid blocks. The correlation between $\log$ span and $\log$ ratio barely changes
($+0.54$ uniform, $+0.48$ at $\rho=32$), so this is a shift of the ratio range, not a change of
structure.

In short: the lower $D$ is a *sampling* effect, as suspected, but the mechanism is the extra freedom
to place the points (which removes a scale bias), not the presence of long arcs. Since every fixed
$\rho$ still has $D_m(\rho)\to D_\infty$, it buys rate of convergence, not a different limit.
Caveat: this diagnostic is at $m=32$; it should be repeated at the largest $m$ we can afford.
