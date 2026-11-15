# Limiting-curve directions from the finite-$m$ shapes

This is an exploratory note, not a new continuum bound or a claim that any of the
models below is correct. The trigger is the $m=100$ restart from the three-piece
IFS candidate: the exact all-pairs discrete distortion starts at $1.98942458399$,
but the optimizer reaches a different part of the shape landscape, at
$1.44282788406$ (unperturbed). This does not refute the IFS on its own; it says the
candidate is not a good guide to the finite-$m$ shapes explored by that optimizer.

## Directions, with tests rather than just names

### 1. Coupled, compatible towers of configurations

Do not choose an independent best local minimum at each $m$ and then infer a curve.
Build $m=2^kN+1$ configurations with nested dyadic parameter values. Optimize at
several levels together, or use a beam of coarse candidates chosen with foresight,
rather than freezing an arbitrary coarse best shape. A useful objective can couple
all levels, for example a weighted maxima of their distortions plus compatibility
penalties. A true limiting curve requires bounded two-sided chord ratios and
Cauchy convergence along nested parameter points.

If vertices are frozen exactly and $D$ remains bounded, the two chord bounds give
a uniform limit and hence a genuine bi-Hölder curve after extending from dyadic
points. This makes nested insertion a construction with a clear success criterion,
not merely a visual interpolation heuristic.

**A preliminary negative test.** Starting with 21 points sampled from the best
stored $m=100$ configuration, inserting midpoints and optimizing **only new
vertices** produced:


| level | $D$ |
|---|---:|
| $m=21$ (initial coarse sampling) | $1.483451068103282$ |
| $m=41$ | $1.4883297046460635$ |
| $m=41 \to 81$ | $1.6086813431453229$ |
| $m=81 \to 161$ | $1.7795343661542535$ |

The old-even-index vertices are preserved exactly, and $D$ was recomputed over all
pairs. This is not the failed canonical two-substep snowflake refinement: the new
vertices were free. But it is one nested annealing run from one coarse choice, not
an optimized joint hierarchy, and not a refutation of all nested constructions.
It says compatibility must influence the selection of coarse shapes. A good test
now is to co-optimize old and new vertices with a limit-compatible norm, or retain
several coarse candidates and select using deeper-level cost.

Files: `scripts/test_nested_insertions.py`, `results/nested_insertion_tower.npz`,
`results/nested_insertion_tower.json`, and `.out`.

### 2. Learn a context-dependent microscopic rule

A single IFS fixes one local pattern globally. The finite-$m$ data instead suggest
isolated stationary configurations with a shared coarse envelope but substantial
small-scale variation. A model with shape states (an arc library), address-dependent
reflections, and possibly slowly changing contraction data has the flexibility to
represent that structure. It could be a finite graph-directed substitution system,
a bounded-turn grammar, or a slowly drifting/non-stationary refinement rule.

The discriminating data are not a few more multistart values of $D_m$. They are:

1. The **restriction-defect spectrum**: the distortion of length-$k$ subarcs of
   larger configurations, compared with the best known $k$-point values.
2. Scale-indexed **turn/phase fields**, measured after realigning coarse geometry.
3. Whether the apparent field is stationary in logarithmic scale, requires context
   from ancestors, or drifts systematically with $m$.
4. Reconstruction tests with the exact all-pairs $D$ and minimum chord separation,
   not only a fitted image.

Large and context-sensitive sub-arc defects would argue against a simple local
substitution grammar and toward a genuinely nonlocal continuum. Small, patterned
defects would support a graph-directed/library model.

### 3. Coarse moduli plus a separate microscopic completion

The stored cloud's dominant deformation has an approximate field
$z(1-z)(a+bz)$. This is a **shape-coordinate model**, not an acceptable final curve:
a smooth polynomial deformation or polygonal completion can be Lipschitz and hence
have infinite snowflake distortion. The failed applications of global Möbius,
cubic-bend, and straight-line deformations likewise do not rule out a coarse field;
they show that applying the field by itself is expensive in $\log D$.

The hypothesis to test is instead
$$\text{finite-scale shape}=\text{coarse profile/mode}+\text{scale-dependent rough completion}.$$
The completion must reproduce measured chord ratios across many $m$, not just small
PCA residuals. A good macro model should become more stable with smoothing and
with continued compatible levels; a micro model should preserve bounded ratios,
not create new near-contacts when extras vertices are inserted.

### 4. Separate the image curve from its index law

The finite configurations record both a point set and a parameter ordering. Some
near-tied shape differences might be substantial in the current coordinates but
nearly disappear after a sensible monotone reparametrization. Test this with a
bounded-slope/DTW alignment, followed by exact chord-ratio evaluation; do not let
an arbitrary alignment create an artificial agreement. If parametric alignment
fails, the shelf has real geometric moduli, not one image with many speeds. A
nonlinear, graph-directed conformal IFS is one way to allow address-dependent macro
geometry without forcing every small piece to be a linearly similar copy of the
whole. Its chord bounds would still have to be checked.

## Important alignment caveat

The `t3.py` decimation-shape comparison in the earlier
`.ultramath/smwh6zmc9/m4/report.md` has an interpolation inconsistency: it uses an
outer grid $0,1,\dots,n-1$ but interpolates on $[0,1]$, rather than interpolating
at the actual index-normalized parameters $j/(n-1)$. This effectively applies an
extra index warp and extrapolation and invalidates its very small shape-agreement
numbers as comparisons in the uniform index parameter. It does **not** affect the
separately recomputed distortions or first-order finite-$m$ certificates. New
cross-$m$ shape claims should use the same normalized index law and no unreported
similarity scaling.

## Accepted scoped structural data: a local-plus-global marginal-contact model

The revised report `.ultramath/smwhdqmht/m8/report.md` was checked and accepted as
**finite-$m$ structural evidence plus conjectural directions**, not a continuum theorem.
In the tested $m=100,200,400$ configurations the distortion computed from adjacent
pairs alone agrees with the full all-pairs distortion to roughly $10^{-14}$. This is a
numerical tie, not a claim that every extremal chord is adjacent. At $m=400$, there
are 82 near-maximal pairs, of which 9 are adjacent and 73 nonadjacent; among 698
near-minimal pairs, 332 are adjacent and 366 nonadjacent. All 301 contiguous
100-vertex blocks of that configuration reproduce the whole distortion to
$1.24\times10^{-13}$. These figures were independently recomputed.

The resulting leading hypothesis is **local-plus-global marginal contact**: the
shortest scale and long-range constraints compete marginally, with neither strictly
dominant. A literal dense $\alpha$-geodesic condition would imply a global planar
isometry of the snowflake interval and is therefore ruled out; the earlier version
of that inference was rejected and removed. The retained continuum interpretation
is a pair of stress measures, one for maximal and one for minimal chord ratios,
with a balance of their log-chord gradients. The proposed next test is weak
convergence of canonically selected first-order multiplier measures along a
shape-compatible branch, taking degeneracy of the multipliers into account.

This does **not** rule out a jointly optimized nested limiting construction. It only
shows that a refinement rule cannot be read directly by decimation of the sampled
branches. The sub-arc-window test was found vacuous as an optimality-deficit test,
which likewise does not eliminate hierarchical models. The global infimum $D_\infty$
does not depend on a selection rule; only a selected branch's limiting distortion,
which may or may not equal $D_\infty$, does.

## Checked compactness and the discrete efficiency gap

The revised report `.ultramath/smwhdqmht/m3/report.md` establishes an elementary
compactness lemma, after two proof-reading corrections. Put $h=1/(m-1)$ and
normalise each configuration using **translated** coordinates
$$w_i=\frac{v_i-v_0}{(m-1)^\alpha r_{\min,m}}.$$
(The translation is essential for the Arzelà–Ascoli boundedness claim; $D$ is
unchanged by it. The node relation is of course in modulus, not a complex number
equal to a real ratio.) At grid points
$$|w_i-w_j|\in[(|i-j|h)^\alpha,D_m(|i-j|h)^\alpha].$$
Uniform-in-$m$ boundedness of the $D_m$ gives boundedness and equicontinuity of
linear-in-index interpolants. A correct coarse estimate from the node bounds is
$$|W_m(x)-W_m(y)|\le \Lambda(|x-y|+h)^\alpha+2\Lambda h^\alpha,$$
with $(|x-y|+h)^\alpha\le|x-y|^\alpha+h^\alpha$; this suffices for equicontinuity.
For fixed $x\ne y$, clipped nearest grid indices are eventually distinct and both
bounds pass to the limit, giving $\gamma$ injective with
$$|\gamma(x)-\gamma(y)|\in[|x-y|^\alpha,\Lambda|x-y|^\alpha].$$
The polyline is only a compactness device; it is not claimed to have finite
continuous distortion.

Restriction of any admissible curve to its uniform nodes gives
$D_m^{\rm free}\le D_m(\gamma)\le D(\gamma)$. Applying the compactness statement
yields $D_m^{\rm free}\to D_\infty$, and every uniform subsequential limit of
near-optimal discrete configurations is optimal. This is conditional on near-
optimality: the stored stationary points/search configurations are not established
to be global $D_m^{\rm free}$ minimisers.

Consequently, for an optimal curve $\gamma$ the discrete free infimum and its own
sample distortions both approach $D_\infty$, so the **necessary** efficiency-gap
condition $D_m(\gamma)/D_m^{\rm free}-1\to0$ holds. It is not sufficient and its
rate is uncontrolled: the three-piece arc's measured gap $0.385$ at $m=100$ does not
exclude it. The only established explicit upper bound in this workspace remains
Koch's $2.096805160875$; the three-piece $2.017935$ is provisional numerical.
The m3 script `kkt_check.py` also reproduces the singleton-pair KKT contradiction
for the non-stationary $m=400/800$ search points (relative residuals $1.0685$ and
$1.0417$).

## Current status

No particular limiting curve is selected. Numerical extrapolation of $D_m$ is not
a reliable substitute for constructing a compatible limiting curve. The independent
run continues on the remaining diagnostics, including turning-angle fields,
parametric-alignment hypotheses and a spectral low-rank shape model. These are
exploratory diagnostics: no continuum certificate is claimed for them.
