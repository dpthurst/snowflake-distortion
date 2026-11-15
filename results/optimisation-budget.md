# What an optimisation to stationarity costs here

Measured on this VM (4 cores), to choose an $m$ for the experiment "converge
one configuration, then rerun the overlap search".

## Cost per SLP iteration

The exact-minimax optimiser is sequential linear programming: each iteration
solves an LP with $2\binom{m}{2}$ rows and $2(m-2)+2$ columns, sparse, about
four nonzeros per row, then does a short line search on the exact objective.
Timed on random problems of that size:

```
   m      LP rows     cols      nnz        LP time
 200       39 800      398     119 000      0.05 s
 400      159 600      798     479 000      0.17 s
 600      359 400    1 198   1 078 000      0.47 s
 800      639 200    1 598   1 917 000      0.94 s
1200    1 438 800    2 398   4 316 000      2.61 s
1600    2 558 400    3 198   7 674 000      4.77 s
```

The line search adds a few evaluations of the dense $\binom{m}{2}$ log-range;
measured on the $\beta$-annealed L-BFGS instead, that is $7$ ms per iteration at
$m=400$, $23$ ms at $800$, $103$ ms at $1600$ — so at $m=1600$ roughly another
second per SLP iteration.

## What that means per configuration

The runs in this workspace took 50–150 SLP iterations per stage and often
needed restarts (the $m=150$ certificate used three outer rounds); budget
$200$–$500$ iterations with restarts:

```
   m      per iteration   one full run      with restarts
 400         0.25 s        1 - 2 min         5 - 15 min
 800         1.2 s         4 - 10 min        20 - 40 min
1600         6 s           20 - 50 min       1.5 - 3 h
```

## Recommendation

**$m=400$ first, then $m=800$ if it pays off.**

* At fixed $m$ the comparison is controlled: the same resolution, before and
  after convergence. That is what isolates "the stored point is not at a local
  optimum" from "finer detail is visible".
* $m=400$ carries the most documented defect: four stored configurations there
  span $D=1.5423$–$1.5526$ ($0.67\%$) and the KKT check fails with relative
  residual $1.07$.
* $m=1600$ costs hours, and its fine structure is inherited from $m=800$ by
  upsampling, so it is the worst value for money as a first test.

A cheaper preliminary at the same $m$: skip the SLP and run the
$\beta$-annealed L-BFGS to a low gradient ($5$–$15$ min at $m=800$). It gives no
stationarity certificate, but if a descent that moves $D$ by a few tenths of a
per cent already changes the overlaps substantially, the point is made.

## Risk

The SLP did not always converge at larger $m$: the earlier lift to $m=800$
ended at `stress=2` instead of a certified point. A failure to certify is a
legitimate outcome of the experiment, not a bug, but it should be budgeted for.
