# Renormalisation for the Koch Hölder constants — an attempt

Goal: prove `rho_min = 2*6^-alpha` and `rho_max = 2*sqrt(211)*69^-alpha`, where
`rho(x,y) = |F(x)-F(y)| / |x-y|^alpha`, `alpha = log3/log4`, for the self-similar
parametrisation `F` of the Koch curve.

**Outcome.** The renormalisation exists, is closed-form, and is exact (Section 3). It
does not by itself give the constants, and Section 5 says why in a precise sense: the
transfer is *critical* (isometric, not contractive), so the constants are a max-plus
eigenvalue, and what is missing is a barrier certificate, not an estimate.

## 1. Setting and the two useful identities

`omega = e^{i pi/3}`, `lambda = 1+omega`, so `lambda*conj(lambda) = 3`, `lambda^2 = 3*omega`,
`lambda/conj(lambda) = omega`. Two-generator IFS, both of ratio `1/sqrt3 = 2^-alpha`:

    phi_1(z) = conj(z)/conj(lambda),      phi_2(z) = 1 - z/lambda,
    F(t/2) = phi_1(F(t)),                 F(1-t/2) = phi_2(F(t)).

Reflection symmetry of the curve (verified to 2e-16, and it is what makes the
orientation-reversing generator tractable):

    F(1-t) = 1 - conj(F(t)).

## 2. Reduced coordinates

Every pair reduces (Lemma 2.1 of the main note) to `x <= 1/2 <= y`. Write `x = u/2`,
`y = 1 - v/2`, `u,v in [0,1]`. Then `|x-y| = delta := 1-(u+v)/2` and

    rho = |E| / delta^alpha,     E = omega*conj(F(u)) + F(v) - lambda.        (2.1)

Both extremal configurations, in these coordinates:

    min: (u,v) = (5/8, 5/8)          delta = 3/8
    max: (u,v) = (63/64, 191/256)    delta = 69/512

## 3. The transfer (exact, all four branches)

State

    Q(a,b,C,d,j; u,v) = | a*F(u) + b*conj(F(v)) - C | / ( d - (u+v)/2^j )^alpha,

so that (2.1) is `Q(omega, 1, lambda, 1, 1; u, v)`. Split `u = (e+u')/2`, `v = (e'+v')/2`
with `e,e' in {0,1}` and relabel `u'' = v'`, `v'' = u'` (harmless: `delta` is symmetric in
`u,v`). Then on the corresponding quarter of the square

    Q(a,b,C,d,j; u,v) = (1/sqrt3) * Q(a'', b'', C'', d'', j+1; u'', v'')        (3.1)

with

    (e,e')   a''      b''      C''                       d''
    (0,0)    b        omega*a  lambda*C                  d
    (0,1)    omega*b  omega*a  lambda*C - b              d - 1/2^(j+1)
    (1,0)    b        a        lambda*C - omega*a        d - 1/2^(j+1)
    (1,1)    omega*b  a        lambda*C - omega*a - b    d - 2/2^(j+1)

The identity is algebraic — it follows from the two IFS relations and the reflection
symmetry — and has been checked numerically along random orbits for depths `n <= 8`:
`max |Q_0 - 3^(-n/2) Q_n| < 3e-14`. In particular `a_n, b_n` forever stay in the six
units `{+-1, +-omega, +-omega^2}`.

Two readings:

* **It is Lemma 2.1 in disguise.** The factor `1/sqrt3` is `2^-alpha` and cancels the
  parameter contraction `1/2`: along a word of length `n`,
  `rho = 3^(-n/2) Q(state_n; u_n,v_n)`. The renormalisation is *neutral*; `rho` is a
  scale-critical quantity, so no contraction is available to produce the constants.
* **The constants are a max-plus eigenvalue.** With `V(s) = sup_{u,v} Q(s;u,v)`, (3.1)
  reads `V(s) = 3^(-1/2) max_{branches} V(s_branch)`. The sharp constants are the
  eigenfunction value of this transfer operator at the initial state.

The normalised state is bounded: numerically `|C_n| / 3^(n/2) -> 1` and `d_n` decreases to
a limit near `1/2`, so in the normalisation `(a,b,C) -> (a,b,C)*3^(-n/2)` the reachable
set is pre-compact. The difficulty is therefore *not* unboundedness.

## 4. The extremal pairs as pre-periodic words

The branch word is the pair of binary digit strings of `(u,v)`:

    min: u = 5/8   = 0.101_2       v = 5/8     = 0.101_2        word (11)(00)(11)(00)(00)...
    max: u = 63/64 = 0.111111_2    v = 191/256 = 0.10111111_2   word (11)(10)(11)(11)(11)(11)(01)(01)(00)...

Both are eventually `(0,0)`, i.e. both extremal pairs are pairs of grid points — the same
statement as Lemma 2.2 of the main note. The minimum sits on the diagonal `u = v` (symmetric
pair; horizontal chord); the maximum does not.

## 5. Why the argument does not close

(i) **Criticality, not unboundedness.** Since the transfer is an isometry of the ratio,
    `V` is a *fixed point* of a max-plus operator, not a limit of an iteration. Proving
    `V <= C_+` therefore needs a **barrier certificate**: an explicit function `Phi >= V`
    with `Phi(s) >= Phi(T s)` for every reachable state `s` and `max Phi <= C_+`. Nothing in
    the self-similarity supplies `Phi`; it has to be built, and the standard way is a
    max-plus linear program on a discretisation of the state space — computer-assisted,
    and not obviously finite.

(ii) **No 1-D reduction.** A tempting simplification is "for each gap `delta` the extremal
    chord is at the symmetric pair `x = 1/2 - delta/2`". This is *false*: at `delta = 1/4`
    the minimal chord is attained at `x = 7/16` with ratio `1/sqrt3`, against `1` at the
    symmetric pair. The problem stays two-dimensional (position and gap); the diagonal is
    extremal only at `delta = 3/8`.

(iii) **The extrema are cusps at grid points.** `R(5/8 +- 3/2000) = 0.4884` against
    `0.4835` at `(5/8,5/8)`. A coarse grid misses the extremum entirely (`N = 41` reads
    `0.5168`), so sampling cannot substitute for the argument — including inside a
    renormalisation scheme, where the extremal words must be identified exactly.

## 6. What would finish it

(a) **Finite extremal-state lemma.** Show that the states which can carry the supremum
    (resp. infimum) of `V` lie in an explicit finite set — conjecturally the pre-periodic
    words of Section 4 together with finitely many neighbours. Then the conjecture reduces
    to finitely many inequalities in `Z[omega]`, which is the illuminating endgame.

(b) **Certified branch-and-bound.** On a box `B` of the reduced square, bound
    `|F(u)-F_N(u)| <= 3^-N` for the level-`N` polygonal approximation and use a Lipschitz
    bound for `Q` on `B` to enclose `rho` rigorously; subdivide until every box away from the
    two extremal points is excluded. Routine, mechanical, and it does not explain anything —
    but it would turn both constants into theorems.

Neither is carried out here.
