# The ladder at the critical exponent $\alpha=1/2$

Whether $([0,1],|\cdot|^{1/2})$ embeds bi-Lipschitz in the plane is the boundary case of the
snowflake question (see the note). If it does, $\sup_m D_m(1/2)<\infty$; if not, $D_m(1/2)$ is
unbounded. `scripts/ladder_half.py` runs the same annealed L-BFGS at $\alpha=1/2$:

```
  m     4        8        16       32       48       64
  D   1.1688   1.4323   1.7987   2.3386   3.4349   4.6863
```

Rising by a factor $\approx1.3$ per doubling, against $1.436\to1.635$ over a factor of 16 at
$\alpha=\log3/\log4$. The last rungs are the least converged (few restarts, long thin
landscape) and this proves nothing either way, but it is what non-existence looks like.
