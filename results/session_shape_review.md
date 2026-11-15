# Verified status of the limiting-shape investigation

This note records what was checked independently, not the unreviewed claims in the
swarm reports. The completed Koch note and Lean file were not changed.

## Mathematical status

* **Near-optimal discretisations have optimal subsequential limits.** The compactness
  argument in `.ultramath/smwhdqmht/m3/report.md` is sound after the explicit
  normalisation/equicontinuity errata noted below. In particular, it is not a claim
  that interpolating polygons themselves have finite continuous distortion. The
  conclusion is conditional on approaching the true finite-$m$ infimum, not on a
  stored local optimizer's value.
* **The proved explicit upper bound remains Koch's $2.096805160875$.** The
  three-piece value $2.0179353894$ is a sampled ratio and remains provisional, as
  the user requested; its global injectivity and two-sided chord bounds were not
  established here.
* **There is no established numerical value or uniqueness of $D_\infty$.** Finite-$m$
  searches constrain $D_m$ in one direction; stored stationary or best-found
  configurations are not automatic global $D_m$ optima. An extrapolated plateau
  above Koch's proved value is inadmissible, but this does not select a lower plateau.
* A compactness clarification used in the working proof is
  $w_i=(v_i-v_0)/((m-1)^\alpha r_{\min,m})$. The origin matters for
  Arzelà–Ascoli boundedness. A sufficient modulus estimate is
  $$|W_m(x)-W_m(y)|\le\Lambda(|x-y|+h)^\alpha+2\Lambda h^\alpha,$$
  rather than the shorter unproved intermediate inequality sometimes printed.

## Finite-$m$ facts that survived checking

### Local and global marginal contact (accepted d7)

At the tested $m=100,200,400$ configurations, adjacent-only chord distortion agrees
with full distortion to about $10^{-14}$. This is a numerical tie, not an assertion
that all extremal pairs are adjacent. At $m=400$, the near-maximal set is
$9$ adjacent $+73$ nonadjacent pairs, and the near-minimal set is
$332$ adjacent $+366$ nonadjacent pairs. All 301 contiguous 100-vertex blocks of
this configuration reproduce the whole $D$ to $1.24\times10^{-13}$.

This supports a **local-plus-global** stress-measure hypothesis, not a dense
$\alpha$-geodesic condition (which would imply an impossible planar isometry). The
value of the proposed next test is canonical first-order multiplier/stress measure
convergence along a shape-compatible branch. Degeneracy of the multipliers has to be
handled explicitly. File: `.ultramath/smwhdqmht/m8/report.md`.

### Cross-ratio and diagonal-gauge diagnostics (accepted d11)

The LP for $D^{red}$ is an unconstrained chord-data relaxation, **not** a necessary
curve-optimality condition. One must not identify it with a realisable planar
reparametrisation. The earlier claim that an arbitrary LP vertex proves high-frequency
slack was retracted: a least-TV optimal solution is much smoother.

The chord cross-ratio fingerprint is exactly invariant under reversal and
conjugation. Full enumeration confirms this, including the signed skewness. It is
therefore not evidence of chirality. The currently measured target band is
$$\text{sd}\approx0.1487\ (m=100),\qquad 0.1604\ (m=400),$$
with skew about $-0.034$ and $-0.026$. These are finite-$m$ classifiers to test for
persistence, not yet necessary continuum conditions. File:
`.ultramath/smwhdqmht/m5/report.md`; `tools/` contains reproduction code.

### Codimension (accepted d13 in restricted scope)

The $\mathbb R^3$ recovery advantage after finite perturbations is real at the tested
budget: for example $m=100$ gave $1.481222849$ in $\mathbb R^3$ versus $1.506574390$
in $\mathbb R^2$, with neither reaching the stored planar optimum. This is not a
bound on global codimension benefit or a proof that planarity is the obstruction.
The next informative test would be a genuinely global $\mathbb R^3$ multistart at
$m=40$, not more claims from local recovery. Finite-distortion plane embedding of
the snowflake interval is already supplied explicitly by Koch; no existence
assumption is needed there. File: `.ultramath/smwhdqmht/m6/report_d13.md`.

### Symmetry and new achieved values (finite-$m$ evidence)

Checked feasible symmetric configurations include
$$D\le1.447066253\ (m=100,\ T_1),\qquad 1.442349749\ (m=100,\ T_2),$$
$$1.497314357\ (m=200,\ T_2),\qquad 1.543697915\ (m=400,\ T_2).$$
They are not stationary certificates. This makes exact symmetry cheap at finite
$m$, but the best $m=100$ configuration is itself the most $T_2$-like of its cloud.
There is no inference from a visually asymmetric cloud to an asymmetric limit.
The claim of a universal minimum over all point-symmetric similarity arcs was not
accepted; even its quoted $D_8$ minimum is contradicted by the full sweep's output
$1.7854661005$ at $w_1=0.96$. This concern is about scope, not the feasibility of
the actual finite-$m$ configurations.

## Two coherent macro-shape branches

Member m1 identified an S-shaped branch, separately from the one-sided arc family.
Independent all-pairs checks and first-order LP calls give

| file | $m$ | recomputed $D$ | multiplier LP |
|---|---:|---:|---|
| `morph400.npz`, `m200 from S(m100)+kick` | 200 | $1.4916292099925$ | not found (one pair on each side) |
| `s_m400.npy` | 400 | $1.5422713964167$ | not found (one pair each) |
| `s_m800_plain.npy` | 800 | $1.5910519607023$ | not found (one pair each) |

Files are under `.ultramath/smwhdqmht/m1/`. These improve on the corresponding
previously stored values but are **not certified stationary and not global optima**.
The failure of the LP is a lack of a certificate here, not a proof of global
non-optimality.

Using endpoint-normalised shapes and the physical index grid (the earlier rejected
interpolation artifact is not used), the distances are approximately
$$d(m_{100}^{\rm old\ best},S_{400})=0.0919,\quad
d(S_{200}^{\rm kick},S_{400})=0.0447,\quad
d(S_{400},S_{800})=0.00652.$$
This is evidence for a reasonably compatible **shape progression**, not exact nested
refinement and not a continuum convergence proof. The old argument that one can infer
a limit only from a single arc-like branch was too strong. The remaining need is
finer S-branch searches at compatible resolutions, followed by first-order/stress
and hierarchical compatibility checks.

## Two methodological warnings for future reading

1. The old interpolation-grid argument comparing selected $m$ shapes was inconsistent:
   it queried $[0,1]$ against sample abscissae $0,1,\dots,m-1$. Its very small
   shape-agreement numbers and its derivative conclusion about sensitivity must not
   be used as though they are uniform-index measurements.
2. The universal leaf-band assertion--$\rho\equiv1$ for every alpha-similarity arc--
   **fails on the saved three-piece grid**. `scripts/check_leaf_band_claim.py`
   evaluates all $39,366$ consecutive distinct-parameter pairs with literal parameter
   differences and obtains
   $$0.9999999995717364\le\rho\le1.0745699320766973.$$
   A maximum example is at parameters $0.4913778731294849$ and
   $0.49137805533128426$. Accordingly it cannot be used as a general structural proof
   that every such arc fails the efficiency-gap test. Different interpolation or
   address discretisations must be explicitly distinguished.

`results/leaf_band_check.out` is the saved output. The rejected three-piece value
stays provisional even if another report calls it a genuine continuum value.

### Additional $m=1600$ S-branch check

The plain file `.ultramath/smwhdqmht/m1/s_m1600_plain.npy` is now available. An
independent all-pairs computation gives $D=1.635530899748628$, with minimum
separation/diameter $0.00199268$. Its imaginary coordinate still has the S-profile,
ranging from $-0.247871$ to $+0.230180$. The endpoint-normalised distances are
$0.004383$ to the plain $m=800$ S configuration and $0.007242$ to $m=400$. These
are further compatible-shape data, not a numerical estimate of $D_\infty$.
The first-order multiplier LP remains infeasible, with one extremal pair on each
side; this is not a stationary-point certificate.

The comparison plot `figures/s-branch-convergence.png` is generated by
`scripts/plot_s_branch_convergence.py`. The joins between sampled nodes are only a
display device, not an admissible polygonal continuum completion.

The perturbed $m=1600$ run and any final m1 report remain to be assessed before a
limiting curve is selected.

### The arc-side comparison changes the interpretation: no unique leading branch

A subsequently completed arc-side $m=1600$ computation produced
$$D_{\mathrm{arc},1600}=1.631583964892305,$$
which is below the plain S-side value $1.635530899748628$ by about $0.24\%$. At
$m=800$, however, the polished arc value $1.593953535872861$ is above the S-side
value $1.5910519607023$ by about $0.18\%$. These are exact all-pairs values on
`.ultramath/smwhdqmht/m1/arc_m1600.npy` and `arc_m800.npy`, not global optima.

The arc is also geometrically compatible with its own earlier configuration:
the endpoint-normalised distance from `arc_m1600` to `arc_m800` is $0.00508758$.
Meanwhile `arc_m1600` is $0.80255$ from the S configuration in the four-variant
comparison. Both arc and S configurations still fail the first-order multiplier LP
(one extremal pair on each side). Consequently, the current evidence supports **two
coherent macro families with a small, resolution-dependent $D$ ordering**, rather
than a unique preferred S branch. The $m=1600$ arc improves the achieved upper
bound for that resolution, but does not certify an optimal branch or distinguish
optimisation from branch selection. The figure generation was updated to use this
newer arc, not the earlier under-polished arc control.

The perturbed S-side value is $1.635914838658449$, with shape distance $0.00345549$
from the plain S-side configuration and no stationary certificate. The kick
## Close-out status of this run

The swarm has been stopped after the model-call budget was exhausted. Four scoped
reports were accepted: marginal stress, compactness/efficiency, diagonal-gauge/CR
diagnostics, and local codimension recovery. The later arc-versus-S configurations
were computed and checked separately from those accepted reports; **no unique
limiting curve or global minimiser has been selected**.

The figure and this note distinguish compatible shapes through several resolutions
from an established continuous limiting curve. The current open tests are matched
first-order/stationary checks and explicit stress/compatibility convergence on each
branch, not another comparison of provisional continuum bounds. The m1 report file
was still a draft when the run closed, with a final arc-side comparison and further
stationary-refinement work assessed independently here. The proved Koch bound
$2.096805160875$ remains the rigorous explicit ceiling; $2.017935$ remains
unverified/provisional.

No further stationary work is claimed as a result. Rejected universal leaf-ratio
assertions, unsupported global exclusions of spectral or hierarchical models, and
invalid use of the sampled three-piece value as a rigorous continuum ceiling are
specifically not endorsed.
