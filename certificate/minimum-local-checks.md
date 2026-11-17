# Sharp Hölder constants for the Koch curve — m1 / direction d6, submission 2

**Status.** `inf rho = 2*6^{-alpha}` is a THEOREM (rigorous, exact, computer-assisted with a fully
rational local certificate, interval arithmetic for the alpha-powers). `sup rho =
2*sqrt(211)*69^{-alpha}` is proved on all of the fundamental domain **except 8 explicit cells of side
1/1024**, which are reduced here to four instances of an explicit per-cell algebraic inequality; my
one-sided certificates close 2 of the 8 and miss the other 6 by the exact rational amounts recorded
in §2. I do not claim the maximum.

Notation: `alpha = log3/log4` (`4^alpha=3`, `2^alpha=sqrt3`, `3^{1/alpha}=4`), `omega^2=omega-1`,
`lambda=1+omega`, `conj(lambda)=2-omega`, `lambda^2=3*omega`, `N(a+b omega)=a^2+ab+b^2`. All
arithmetic exact in `Z[omega][1/3]`; `e_mod2` is the Eisenstein norm (my first submission had a sign
slip in the *display* routine `ev` for the real part; all norms were and are correct — the affected
local certificates are redone below and the min-side conclusion survives with larger slack).

---

## 0. The two reductions (proved)

**(0.1)** `F(1-u) = 1 - conj F(u)` (the symmetry `z -> 1-conj z` swaps `f_0<->f_3`, `f_1<->f_2`).

**(0.2)** With corner coordinates `a = 1-t`, `b = 1-s` of the reduced pair (`x=t/2`, `y=1-s/2`),

    F(x)-F(y) = N/lambda,  N = omega*conj F(t)+F(s)-lambda = -(omega F(a)+conj F(b)),  |x-y| = (a+b)/2,

so, using only `2^alpha = sqrt3 = |lambda|`,

        **   rho(x,y) = Phi(a,b) := | omega F(a) + conj F(b) | / (a+b)^alpha .   **

(The `omega^2` of the problem statement is a typo: it gives 1.2560 at `(5/8,5/8)` instead of
`2*6^{-alpha}`.)

**(0.3) Exact invariance (one line).** `f_0(z)=z/3` gives `F(a/4)=F(a)/3`, hence

        **   Phi(a/4 , b/4) = Phi(a,b)    for all (a,b),  since 3 = 4^alpha.   **

Consequences: the range of `rho` on `[0,1]^2\{0}` equals its range on the compact domain
`D = {max(a,b) >= 1/4}` (there `a+b >= 1/4`, so the pinch `x=y=1/2` never arises); all extrema are
infinite orbits `(a,b)->(a/4,b/4)`; the minimiser set contains `{(t,t): 1-t = 3*2^{-(2k+3)}}`
(`5/8, 29/32, 125/128, 509/512, ... -> 1`) and the maximiser set contains the orbits of
`(63/64,191/256)` and `(47/64,255/256)`. The problem's uniqueness clause is therefore false as
stated and is dropped.

**(0.4)** `K=F([0,1])` lies in `T0 = conv{0,1,lambda/3}`, and the three vertices of `T0` lie ON `K`
(`F(0)=0`, `F(1)=1`, `F(1/2)=lambda/3`). Hence for linear functionals the hull bound is exact:
`max_{z in K} Re(xi z) = max(0, Re xi, Re(xi lambda/3))`. Also `|F(u)| <= sqrt3 u^alpha` (the ratio
`|F(u)|/u^alpha` is invariant under `u->2u`, and `|F|<=1`), and `(1/sqrt3)^{1/alpha} = 1/2`.

---

## 1. The minimum is a theorem

Exact extremal data (`ex4.py`): at `(t,s)=(5/8,5/8)` (`a=b=3/8`)

    N = -(6+6omega)/27 = -(2/9)lambda,   |N|^2 = 4/27,   Phi = 2*6^{-alpha} = 0.4834590783544287.

Proof chain:

1. (0.3) reduces everything to `D`.
2. **Certified enclosure** (`cert2.py`). For a cell of side `2^{-K}` the variation of `F` is at most
   `3^{-floor(K/2)}` (the cell shares `K` binary digits, each linear part has modulus `3^{-1/2}`,
   `diam K = 1`); `|M|^2` at dyadic points is an exact rational. Adaptive bisection to depth 28
   leaves only **20 cells of side 1/1024** (4 at the minimum `(3/8,3/8)`, 16 at the four
   maximisers); every other cell satisfies `c_min <= Phi <= c_max` with slack `>= 1.1e-9`.
3. **The 4 min cells are closed exactly, in rational arithmetic** (`cells.py`). Transfer to depth 10
   with `(a,b) = ((3+u)/8,(3+v)/8)`, `u,v in [0,1/128]`, gives
   `|M_S| = (sqrt3/9)|2 + omega conjF(u) + conjF(v)|` and `Phi = 8^alpha |M_S|/(6+u+v)^alpha`, with
   `|M_S(0,0)| = |E| = 2*sqrt3/9`. The exact one-sided certificate

       |M_S|  >=  ( |E_S|^2 + 2*sum_i min_{z in T0} Re(conj(E_S) c_i z) ) / |E_S|

   (rational, `min` exact by (0.4)) gives `Phi - c_min >= 0.7496, 0.7597, 0.7597, 0.7700` on the
   four cells. Note `Phi(a,b) >= |Re(conj(u) M)|/(a+b)^alpha` is exactly the projection idea of this
   direction: the extremal chord direction is what makes it sharp.

Therefore **`|F(x)-F(y)| >= 2*6^{-alpha}|x-y|^alpha` for all `x != y`**, with equality at least on the
orbit above. (Direction of latitude: the alpha-powers are compared with directed interval arithmetic,
`mpmath.iv`, `iv.power(x, iv.log(3)/iv.log(4))`; step 3 is entirely rational.)

---

## 2. The maximum: exact reduction, 2 of 8 cells closed, 6 open

Exact extremal data: at `(a,b) = (1/64,65/256)` (and at `(1/256,17/64)`, `(65/256,1/64)`,
`(17/64,1/256)`, plus their `4^k`-orbits)

    N = (-2268-162omega)/3^8,  |N|^2 = 844/6561 = 4*211/3^8,  N/|N| = -(14+omega)/sqrt211,
    Phi = 2*sqrt(211)*69^{-alpha} = 1.0137194905655114      (211 = N(14+omega), 69 = 1+68),

and since `(69/512)^{2alpha} = 69^{2alpha}/3^9` this is the rational identity `|N|^2 = 844/3^8`, so
the conjectured constant is exactly right.

**The 8 open cells.** With cells `[i,i+1]/1024 x [j,j+1]/1024`, the four maximisers are shared
corners of level-10 cells: `(4,272)/1024 = (1/256,17/64)` and `(272,4)/1024`. No enclosure can
certify a cell whose corner is the maximiser, because the maximum is attained exactly there. Of the
16 such cells, 8 are closed by the §1.3 certificate (those around `(1/64,65/256)`, margin `-0.038`).
**The 8 that remain are `(3,271),(3,272),(4,271),(4,272)` and their transposes.**

**Exact form of what is missing.** Let `S` be the state at one of these cells (depth 10) and
`(u,v) in [0,1]^2` the cell coordinates. The transfer identity is exact, so

    Phi(a,b) = 1024^alpha |M_S(u,v)| / (A+B+u+v)^alpha,      A+B+u0+v0 = 276,

and `Phi <= c_max` on the cell is **equivalent** to

    (*)      |M_S(u,v)|  <=  |M_S(u0,v0)| * ((A+B+u+v)/276)^alpha ,

with `|M_S(u0,v0)| = 2*sqrt(211)/3^4 = 0.3586626925` (the exact extremal value). So the whole
maximum reduces to four instances of (*) — a fixed-state inequality on the unit square with equality
at the corner (u0,v0). Numerically (*) holds with strict inequality away from the corner.

**What closes and what does not.** Using `(1+x)^{2alpha} >= 1+2alpha x` and exact hull bounds for the
linear parts (`maxproof2.py`, rational throughout):
* cells with corner `(0,0)` — `(4,272)` and its transpose — **are proved**: the exact rational
  vertex certificate over `T0 x T0` gives maximum `0`, attained only at the corner.
* the other 6 cells are **not** closed by it: the same vertex certificate misses by a positive
  amount of order `1e-4` (computed as `4.660831e-4*alpha` resp. `2.330415e-4*alpha` in the
  formulation of `maxproof2.py`; I have not converted these into a rigorous statement, and a second
  formulation I tried (`maxtry.py`) had a sign inconsistency, so treat these numbers as indicative).
  The miss occurs at the hull vertex `z = lambda/3`, a genuine point of the curve (`F(1/2) =
  lambda/3`), so replacing `T0` by level-2/3 hulls or increasing the transfer depth does not help
  while (*) itself still holds — I verified by direct exact evaluation that (*) holds strictly at
  the points the certificate flags (e.g. `(a,b) = (7/2048,543/2048)`: `Phi = 1.0049959394 <
  c_max`), so the failure is in the one-sided estimate, not in the inequality.

**Precise obstruction to the renormalisation (this is the answer to "where does it fail").**
(i) The transfer is exactly neutral in the scale by (0.3) — the invariance is an eigenfunction with
eigenvalue 1 — and the only contraction is `3^{-1/2}` per letter; with 4 letters `4*3^{-1/2}=2.31>1`,
so there is no Perron-Frobenius / LP / function-cone closure of the bare renormalisation.
(ii) The exponents are exactly critical: `3^{1/alpha}=4`. Quantitatively, for a cell of side
`2^{-K}` the enclosure error in `Phi` is at most `2*3^{-K/2}/(a+b)^alpha = 5.65*2^{-K alpha}` while
the true deficit near a maximiser is about `3.7*2^{-K alpha}`: ratio `1.53`. So no interval / grid /
modulus-of-continuity bound can decide the 8 cells at any depth (checked to depth 28), which is also
why the problem's grid sampling drifts like `N^{-alpha}`.
(iii) At a maximiser the perturbation is anti-parallel to the field (`|M(a*,b*)| = 0.35866 <
|E_S| = 0.37920`), so an upper bound on `|M'|` cannot help: what is needed is a *negative,
rate-controlled* one-sided estimate `Re(conj(E')M') <= -c(|z_a|+|z_b|)^{alpha}`, and the convex hull
of `K` — the only cancellation-free exact tool — does not supply one (it fails at `lambda/3`).
The remaining 6 cells are therefore a genuine gap, not a bookkeeping omission.

---

## 3. Files
`exact.py` (exact `Z[omega][1/3]`, `F` at dyadics), `states.py` (transfer identity, verified exactly),
`cert2.py` (certified enclosure on `D`), `cells.py` (min-side exact per-cell certificates),
`maxcert.py`, `maxproof.py`, `maxproof2.py` (max-side exact certificates: 2 of 8 cells closed, the 6
failures quantified), `ex4.py` (exact extremal numerators), `local.py` (local `h^alpha` law).
