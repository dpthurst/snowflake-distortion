# 100-point restart from the three-piece candidate

The restart samples the **parameter interval**, not the mixed-radix IFS grid indices:
$$t_j=j/99,\qquad v_j=\gamma(t_j),\quad j=0,\dots,99.$$
The coordinates $\gamma(t_j)$ are evaluated directly from the saved IFS coefficients in
`.ultramath/smwh6zmc9/m4/coords_3piece_sym.npz`. All 100 are distinct.

The discrete objective is
$$D(v)=\max_{i<j}r_{ij}/\min_{i<j}r_{ij},\qquad
r_{ij}=|v_i-v_j|/|i-j|^\alpha,\qquad \alpha=\log3/\log4.$$
The unperturbed start has $D=1.989424583988947$ (all pairs recomputed). It is not the
old erroneous value $18.55$, which arose from an irregular parameter sample and an
incompatible fixed-vertex gauge. Here the optimizer input is first regauged so its
anchor $v_1=1$ is consistent; reporting and shape comparisons use $v_0=0$, $v_{99}=1$.

## Results

The algorithm is a $\beta=0.1,0.3,1,\dots,10^6$ annealed L-BFGS run, three passes at
each $\beta$, followed by a guarded sequential-LP polish. The relative perturbation is
$10^{-5}$ times the median initial step length, independently normal for real and
imaginary parts of non-endpoint vertices.

| start | final $D$ | tight maximum pairs | tight minimum pairs | multiplier residual (sup norm) |
|---|---:|---:|---:|---:|
| unperturbed | **1.4428278840594533** | 24 | 167 | $1.03\times10^{-14}$ |
| $10^{-5}$ perturbation, seed 41 | **1.4439157038101809** | 28 | 164 | $1.09\times10^{-13}$ |
| $10^{-5}$ perturbation, seed 73 | **1.4438698021948315** | 28 | 168 | $2.91\times10^{-15}$ |

A separate explicit Python loop over all $4950$ pairs reproduced the first value as
$1.442827884059453$ and the other two to the same precision. Each final shape has 100
distinct vertices. The last run's stopping point was verified independently using
the non-smooth first-order multiplier linear program at tolerance $10^{-8}$ in log-ratio.
These are numerical local stationarity checks, **not global optimality certificates**.

## Interpretation

The restart lands on the familiar low-$1.44$ shelf, but not at the best known $m=100$
configuration, whose value is $1.4362283999316854$. The unperturbed result is about
$0.46\%$ higher in distortion. The normalized sup-distances to that best shape are
$0.52245$, $0.52915$, and $0.53788$; they are not the same shape up to the four usual
variants (conjugation and/or parameter reversal).

For context, all three are nearest, among the 81 stored configurations in
`.scratch/final_m100.npz`, to stored configuration 9, whose old value is $1.4392148630$.
Their normalized distances to configuration 9 are $0.03050$, $0.04561$, $0.05761$.
Thus this is not a return to the starting IFS shape, nor an exact recovery of the old
best point. It is an $m=100$ discrete optimization result, not a new bound for the
continuum candidate.

## Reproducibility

Run from the project root:

```bash
python3 -u scripts/restart_new_candidate_m100.py
python3 -u scripts/polish_new_candidate_restart_m100.py
```

The second command continues the seed-41 polish after its initial 400 steps, since
that run had not yet reached a first-order multiplier certificate.

- `results/new_candidate_uniform_m100_restarts.npz`: sample, annealed and polished
  coordinates, multiplier arrays.
- `results/new_candidate_uniform_m100_restarts.json`: numerical metadata, tight pairs,
  comparison values and independent recomputations.
- `results/new_candidate_uniform_m100_restarts.out`: annealing/polish log.
- `results/new_candidate_uniform_m100_additional_polish.out`: extra seed-41 polish.
- `results/new_candidate_uniform_m100_independent_checks.out`: separate all-pairs and
  stationarity checks.
