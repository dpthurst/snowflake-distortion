# Is any finite-$m$ distortion known to be a global minimum?

**No.** Every value $D_m$ quoted in this workspace is the best configuration
found, with at most a first-order stationarity certificate. Nothing here proves
that a better configuration does not exist at any $m\ge4$.

## What is certified, and what that certificate is

* $m=400$ and $m=800$ (today): first-order KKT points of the minimax — a
  multiplier vector on the tight max and min pairs exists (`stress` status 0),
  with $D=1.5417920000$ and $1.5900261898$. The stored points from which they
  were obtained failed the same test.
* $m=40$, $m=100$ (earlier sessions): exact stationary points, many of them —
  $42$ distinct stationary values at $m=100$ from the census runs.
* Stationarity says no descent direction exists. It says nothing about other
  basins.

## What the multistart data actually show

Clustering the stored candidates at $1\times10^{-4}$:

```
        candidates   distinct clusters   best        runner-up   gap      spread
m=40        51              31            1.352310    1.353037    7.3e-4   0.0164
m=100       81              40            1.436228    1.436991    7.6e-4   0.0129
```

So the landscape at small $m$ is crowded: dozens of distinct local minima within
about one per cent of each other, with the best isolated by only $0.05\%$. The
best cluster is also small (one or two candidates), so it is not a large basin
that many starts fall into — repeated agreement on the same best value is
evidence about how often that basin is *found*, not evidence that no better one
exists.

## Why there is no floor

A global claim needs a rigorous lower bound: a proof that every configuration of
$m$ points has $D$ at least some value. Nothing in this workspace provides one.
The only rigorous statements are

* $D_m\ge1$ trivially, and $D_3=1$ exactly: three points at indices $0,1,2$ can
  be placed so that all three ratios equal one (a triangle with sides
  $c,c,2^{\alpha}c$ exists because $2^{\alpha}\le2$). So the first nontrivial
  case is $m=4$.
* Koch's $2.096805160875$ is an *upper* bound for the continuum, not a lower one.

The "floors" direction tried earlier (minimise the distortion of an index subset
over planar configurations, using Cayley–Menger) is rigorous in principle but
every numerical evaluation of it returns an achieved configuration, i.e. an
upper bound on the floor and hence no lower bound on $D_m$ at all without
interval arithmetic.

## Best found at small $m$ (40 starts each, not certificates)

```
m=4   1.062766196
m=5   1.093091572
m=6   1.128130050
```

## What could actually be certified

$m=4$ is four-dimensional after the similarity gauge ($v_0=0$, $v_3=1$ leaves
two free complex points). A rigorous interval branch-and-bound on a box in
$\mathbb{R}^4$ — discard a box when $\underline{\max r}/\overline{\min r}$
exceeds the best known value, subdivide otherwise — would settle $m=4$, and
probably $m=5$ (six dimensions) on this machine. Beyond that the dimension
$2(m-2)$ grows too fast to subdivide.
