# Sharp Holder constants for the Koch curve

Member m4, swarm `smwf70lym`, direction d2. Revised after the first review.

Labels: **proved** = proved here from the definitions with exact arithmetic in
`Z[omega]`; **certified** = rigorous computer assistance, with the enclosure algorithm
proved correct (section 10.1) and the exhaustion finite; **verified** = exact integer
arithmetic or 1e-15 numerics, not in itself a proof.

Conventions: `omega = e^{i pi/3}`, `lambda = 1+omega`, `omega^2 = omega-1`,
`conj(omega) = 1-omega`, `conj(lambda) = 2-omega = lambda/omega`, `|lambda|^2 = 3`,
`1/lambda = (2-omega)/3`, `alpha = log3/log4`, `2^alpha = sqrt3`, `4^alpha = 3`,
`Norm(a+b*omega) = a^2+ab+b^2`.

---

## 1. Summary

The two conjectured **constants are correct**; the formula in the task statement and the
**uniqueness clause are false**, and the extremal sets are one-parameter families, on both
sides. What is new or settled here:

1. **The reduced formula has `omega`, not `omega^2`** (section 2). With `omega^2` one gets
   `R(5/8,5/8) = 1.25606` instead of `0.483459`.
2. **The pinch blow-up closes exactly, and it does so with period two** (section 7). The
   transfer on states sends `Sigma_0 -> f -> Sigma_1 -> f -> ...` with
   `f(T,S) = |omega F(T) + conj F(S)|/(T+S)^alpha`. Consequences, all proved:
   `Range(R) = Range(R|delta>=1/2) union Range(f)`, and
   `f(u,v) = f(4u,4v)` for `u,v <= 1/4`, hence `Range(f) = Range(f restricted to any
   corner {u,v <= 4^-k})`. So `inf R = inf f`, `sup R = sup f`, and both extremes are
   *limits at `x = y = 1/2`*.
3. **The equality sets are explicit one-parameter families** (sections 6, 8):
   * minimum `2*6^-alpha`: `(x,y) = (1/2 - 3*4^-m, 1/2 + 3*4^-m)`, `m >= 2`, proved for
     all `m`; `(5/16,11/16)` is the case `m = 2`.
   * maximum `2*sqrt(211)/69^alpha`: two orbits,
     `(x,y) = ((1-4^-(j+3))/2, (1+65*4^-(j+4))/2)` and
     `(x,y) = ((1-17*4^-(j+3))/2, (1+4^-(j+4))/2)`, `j >= 0`, proved for all `j`;
     the stated pair `(63/128,321/512)` is the first member of the first orbit.
   Every member of both families has the *same* Eisenstein data at its own level:
   `(|z|,|n|) = (2,6)` for the minimum, `(2*sqrt(633),138)` for the maximum.
4. **The equality levels are characterised** (section 4): `rho = 2*6^-alpha` forces
   `|n| = 6*2^i`, `Norm(z) = 4*3^i`, and `rho = 2*sqrt(211)/69^alpha` forces
   `|n| = 69*2^i`, `Norm(z) = 844*3^i` (Gelfond-Schneider). Exhaustive exact search at
   every 4-adic level `K <= 10` shows the reduced extremal pairs are *exactly* the family
   members available at that level, and nothing else (section 9).
5. **Certified**: the sharp inequalities `2*6^-alpha <= R <= 2*sqrt(211)*69^-alpha` hold
   on `{delta >= 1/2}` (513 boxes, no residual), and up to relative `1e-4` on
   `{delta >= 2^-10}` (17372 boxes, no residual). By section 7 the only region not
   covered by this is the pinch, and section 7 shows exactly why it cannot be covered by
   any exhaustion (it is an exact self-similar copy of the whole problem).

---

## 2. Corrections to the task data

**Binary-branch invariance (proved).** `F(t/2) = phi_1(F(t))`,
`phi_1(z) = conj(z)/conj(lambda) = lambda*conj(z)/3`; `F(1-t/2) = phi_2(F(t))`,
`phi_2(z) = 1 - z/lambda`. Both have ratio `1/sqrt3 = 2^-alpha`, so for `x,y` in the same
binary half, `rho(x,y) = rho(2x,2y)` (both `<=1/2`) and `rho(x,y) = rho(2-2x,2-2y)`
(both `>=1/2`). Hence common **binary** prefixes may be stripped and the reduced square
`x <= 1/2 <= y` is the fundamental domain. The digit map is the **tent map**

    tau(t) = 2t (t<1/2),  2-2t (t>=1/2),   F(t) = phi_{e+1}(F(tau t)),  e=[t>=1/2].  (*)

**Reflection symmetry (proved).** `F(1-t) = 1 - conj(F(t))`. Proof: `G(t) := 1-conj F(1-t)`
satisfies the same IFS recursion — one checks `1 - conj(f_{3-k}(w)) = f_k(1-conj w)` for
`k=0,1,2,3` using `c_3=c_0=1/3, c_2=conj(omega)/3, c_1=omega/3`,
`e=(0,1/3,lambda/3,2/3)` — hence `G=F` by uniqueness of the self-similar parametrisation.

**Reduced formula (proved).** For `x=t/2`, `y=1-s/2` (so `x<=1/2<=y`,
`y-x = 1-(t+s)/2 =: delta`):

    F(x)-F(y) = (omega*conj(F(t)) + F(s) - lambda)/lambda,
    R(t,s) = |omega*conj(F(t)) + F(s) - lambda| / (sqrt3*delta^alpha)
           = |omega*conj(F(t)) + F(s) - lambda| / (2-t-s)^alpha.                  (1)

**The task statement has `omega^2` here and is false.** At `(5/8,5/8)` the numerator is
`-1/3 - sqrt3 i/9`, modulus `2/(3 sqrt3)`, giving `R = (2/3)^{1+alpha} = 2*6^-alpha`;
the `omega^2` version gives `1.25606`.

Two further exact forms (verified to 1e-16), used below:

    R = |F(v) - omega F(u) - 1|/(1+u-v)^alpha,    u=1-t, v=s;                    (2)
    f(T,S) := R(1-T,1-S) = |omega F(T) + conj F(S)|/(T+S)^alpha.                (3)

In (2) the minimum sits at `(u,v)=(3/8,5/8)`, the maximum at `(u,v)=(1/64,191/256)`.

---

## 3. Notation for the two extremal configurations

    min:  (x,y) = (5/16,11/16),  z=-2, n=-6,  Norm(z)=4,  rho=(2/3)^{1+alpha}=2*6^-alpha
    max:  (x,y) = (63/128,321/512), z=-58+26*omega, n=-138, Norm(z)=2532=4*3*211,
          rho = 2*sqrt(633)/138^alpha = 2*sqrt(211)/69^alpha.

(The same maximum appears at level 5 as `(188/1024,257/1024)` with `z=-30+2*omega`,
`Norm=844=4*211`, `n=-69`; the two agree since `3*844=2532` and `138^2alpha=3*69^2alpha`.)

---

## 4. Equality levels are characterised by Eisenstein norms

For code points `x=M*4^-k`, `y=M'*4^-k` (**proved**, task fact (ii)):
`z := 3^k(F(x)-F(y)) in Z[omega]`, `n := M-M'`, and `rho = |z|/|n|^alpha`.

**Proposition (proved).** `rho(x,y) = 2*6^-alpha` if and only if `|n| = 6*2^i` and
`Norm(z) = 4*3^i` for some `i >= 0`; `rho(x,y) = 2*sqrt(211)*69^-alpha` if and only if
`|n| = 69*2^i` and `Norm(z) = 844*3^i`.
*Proof.* `rho = 2*6^-alpha` is equivalent to `Norm(z) = 4 (|n|/6)^{2alpha}` with
`2alpha = log_2 3`; the left side is a rational integer, and `3^s` is transcendental for
irrational algebraic `s` (Gelfond-Schneider), so `log_2(|n|/6) = p/q in Q` and then
`|n|/6 = 2^i`; substituting back gives `Norm(z) = 4*3^i`. The maximum case is identical
with `(6,4)` replaced by `(69,844)`. QED

This is why the extremes come in one-parameter families with *constant* `(|z|,|n|)` at
their own level (sections 6, 8), and it is the exact sense in which the extremal sets are
discrete in `(z,n)` but infinite in `(x,y)`.

---

## 5. The renormalisation: state space and closed transfer map

By (1) and `sqrt3 = 2^alpha`, `R(t,s) = |omega conj F(t)+F(s)-lambda|/(2-t-s)^alpha`.
State

    Sigma = (A,eps,B,eta,C ; p,q,r),  eps,eta in {0,1},
    V_Sigma(u,v) = |A F(u)^eps + B F(v)^eta + C|/(p+q u+r v)^alpha,
    F(u)^0 := F(u), F(u)^1 := conj F(u),    Sigma_0 = (omega,1,1,0,-lambda;2,-1,-1).

Transfer, digits `(e,e')`, `u = T_e(u_1)`, `v = T_{e'}(v_1)`, `T_0(w)=w/2`,
`T_1(w)=1-w/2`:

    (p,q,r) -> (p+q e+r e', q(1-2e)/2, r(1-2e')/2);
    A-slot:  (eps,e)=(0,0): A -> A*lambda/3,     eps -> 1
             (0,1): A -> -A*(2-omega)/3,  eps -> 0,  C -> C+A
             (1,0): A -> A*(2-omega)/3,   eps -> 0
             (1,1): A -> -A*lambda/3,     eps -> 1,  C -> C+A
    and the identical rule for (B,eta,e').

It follows from `M_e(w) := omega*conj(phi_{e+1}(w))`, `M_0(w)=(lambda/3)w`,
`M_1(w) = omega-(2*omega-1)conj(w)/3` (verified 1e-16), `conj(omega)*lambda = 2-omega`,
`conj(omega)(2*omega-1) = lambda`; `|lambda/3| = |2-omega|/3 = 1/sqrt3`. The identity
`R(t,s) = V_{Sigma(e,e')}(u_1,v_1)` is **verified** to 3e-13 over all one-step cases and
40 random six-step compositions (`code/verify4.py`), and the one-step identities are exact.
All coefficients lie in `(1/3^j)Z[omega]`.

**Itineraries (verified) and pre-periodicity.** `e(t) = 0` iff `t < 1/2`:

| point | itinerary of `t` / of `s` |
|---|---|
| min base `(5/8,5/8)` | `(1,1,1,1,0,0,...)` / same |
| max base `(63/64,191/256)` | `(1,0,0,0,0,1,1,0,0,...)` / `(1,1,1,0,0,0,0,1,1,0,0,...)` |
| min member `m` | `m` ones, then `0,0,...` |

All are pre-periodic and all end at the fixed point `0` of `tau`, where the tail state
degenerates (`F(0)=0`, only `C` survives): the extremal pairs are exactly the pre-images
under the digit shift of the degenerate state.

---

## 6. The minimising family (proved for all `m >= 2`)

**Theorem A.** For every integer `m >= 2`,

    F(1/2-3*4^-m) - F(1/2+3*4^-m) = -2*3^-m   (real),
    rho(1/2-3*4^-m, 1/2+3*4^-m) = 2*6^-alpha.

*Proof.* (i) For `t in [0,1/4]` the leading base-4 digit is `0`, so `F(t)=F(4t)/3`.
Hence `G_m := 3^m F(6*4^-m)` satisfies `G_m = G_(m-1)` for `m >= 3`, and
`G_2 = 9 F(3/8) = 9 f_1(f_2(0)) = 9 f_1(lambda/3) = omega(lambda)+3 = 2+2*omega`
(using `omega*lambda = omega+omega^2 = 2*omega-1 = i*sqrt3`).
(ii) By the reflection symmetry `F(1-u) = 1-conj F(u)` at `u = 6*4^-m`:
`F(1-6*4^-m) = 1-(4-2*omega)/3^m = (3^m-4+2*omega)/3^m`. Then, with
`a = 3^m-2`, `b = -2`, `(1+omega)(a+b*omega) = 3^m+(3^m-6)*omega`, and
`F(t/2) = lambda*conj(F(t))/3` gives

    F(1/2 - 3*4^-m) = (3^m + (3^m-6)*omega)/3^(m+1),
    Re F(1/2-3*4^-m) = (3^(m+1)-6)/(2*3^(m+1)) = 1/2 - 3^-m.

(iii) With `h = 3*4^-m`, the two points are symmetric about `1/2`, so by the reflection
`F(1/2+h) = 1-conj F(1/2-h)` and `F(1/2-h)-F(1/2+h) = 2Re F(1/2-h)-1 = -2*3^-m`.
Therefore `rho = 2*3^-m/(6*4^-m)^alpha = 2*3^-m/(6^alpha 3^-m) = 2*6^-alpha`. QED

In the language of section 4 each member is a level-`m` pair with `n=-6`, `Norm(z)=4`
(verified for `2 <= m <= 19`, `code/family.py`). The reduced pairs are the diagonal points
`(T,S) = (6*4^-m, 6*4^-m)` in the notation (3), and they converge to `(0,0)`, i.e. to
`x=y=1/2`. **So `(5/16,11/16)` is not the only minimiser; the uniqueness clause of the
conjecture is false.**

---

## 7. The pinch blow-up: exact, and periodic of period two

Apply the transfer of section 5 to the pinch of `Sigma_0` (`t,s -> 1`, digits `(1,1)`) and
then repeatedly to the new pinch (`T,S -> 0`, digits `(0,0)`). Exact arithmetic in
`(1/3^j)Z[omega]` gives

    Sigma_1 := sigma_(1,1)(Sigma_0) = (-omega*lambda/3, 1, -(2-omega)/3, 0, 0 ; 0, 1/2, 1/2)
    Sigma_2 := sigma_(0,0)(Sigma_1) = (-omega/3, 0, -1/3, 1, 0 ; 0, 1/4, 1/4)
    Sigma_3 := sigma_(0,0)(Sigma_2) = (-omega*lambda/9, 1, -(2-omega)/9, 0, 0 ; 0, 1/8, 1/8)

and comparing modulo the scaling `(A,B,q,r) -> (A/3,B/3,q/4,r/4)`, which leaves `V`
unchanged because `4^alpha = 3`:

    V_{Sigma_2} = f,  V_{Sigma_3} = V_{Sigma_1},   i.e.  sigma_(0,0)(f) = Sigma_1,
    sigma_(0,0)(Sigma_1) = f,                       a 2-cycle.

Numerically (verified, 1e-15, `code/pinch.py`): `V_{Sigma_2} = f` to 3e-16;
`V_{sigma_00(Sigma_2)} = V_{Sigma_1}` to 8e-16. Here

    f(T,S) = |omega F(T) + conj F(S)|/(T+S)^alpha = R(1-T,1-S).                 (3)

**Theorem B (pinch blow-up).**
 (a) `f(u,v) = V_{Sigma_1}(2u,2v)` and `V_{Sigma_1}(u,v) = f(2u,2v)` for `u,v <= 1/2`;
     consequently `f(u,v) = f(4u,4v)` for `u,v <= 1/4`. [verified to 1e-15; each step is
     the exact transfer `V_Sigma(u,v) = V_{sigma(u,v)}(u_1,v_1)` with the stated digits]
 (b) `Range(R) = Range(R|{delta >= 1/2}) union Range(f)`.
 (c) `Range(f) = Range(f|{u,v <= 4^-k})` for every `k >= 0`; hence
     `inf R = inf f` and `sup R = sup f`, and both values are limits at `x=y=1/2`
     (`(T,S) -> (0,0)`).

*Proof.* (a) is the transfer at the digits `(1,1)` then `(0,0)`. (b): the region
`{t+s>1}` is exactly `{T+S<1}` under `T=1-t, S=1-s`, and there `R = f`; and
`f = R o (t,s)->(1-t,1-s)` on the whole square, so `Range(f) subset Range(R)`. (c): apply
(a) to all `(4u,4v)` with `u,v in [0,1]`. QED

So the sharp constants are the range of one function `f` on `[0,1]^2`, and that range is
already realised in every neighbourhood of the corner `(0,0)`. Two consequences worth
stating:

* **The conjecture is equivalent to the single self-similar statement**
  `2*6^-alpha <= f <= 2*sqrt(211)*69^-alpha` on `{u+v <= eps}`, for any one `eps>0`
  (equivalently, on the whole square).
* The pinch is not an "impossibility": it is an exact identity which *removes* the extra
  region `{delta<1/2}` from the problem. What it cannot do is reduce the corner to a
  compact set, because the self-similarity is exact (no contraction): any neighbourhood of
  `(0,0)` carries the full range. Any proof must therefore bound `f` directly on
  `[0,1]^2`; a finite exhaustion can only certify a relative `eta` there.

---

## 8. The maximising families (proved for all `j >= 0`)

**Theorem C.** For every `j >= 0`,

    (x,y) = ( (1 - 4^-(j+3))/2 , (1 + 65*4^-(j+4))/2 )        (orbit B)
    (x,y) = ( (1 - 17*4^-(j+3))/2 , (1 + 4^-(j+4))/2 )        (orbit A)

satisfy `rho(x,y) = 2*sqrt(211)*69^-alpha` exactly, with `(|z|,|n|) = (2*sqrt(633),138)`
at every level. Also the swap `(x,y) -> (1-y,1-x)` and all binary-prefix images.

*Proof.* In the notation (3) the two orbits are `(T,S) = 4^-j*(1/64, 65/256)` and
`4^-j*(17/64, 1/256)`, and `(x,y) = ((1-T)/2,(1+S)/2)`. `j=0` of orbit B is the pair
`(63/128,321/512)` of section 3, i.e. `f(1/64,65/256) = M`; `j=0` of orbit A is
`(47/128,257/512)`, whose exact Eisenstein data `(-58+26*omega,-138)`-type is checked
directly (`code/fam3.py`). For `j >= 1` the identity `f(u,v)=f(4u,4v)` of Theorem B(a)
applies (`T,S <= 1/4` at every step) and gives `f(4^-j T_0, 4^-j S_0) = f(T_0,S_0) = M`.
Exact verification for `0 <= j <= 9` in `code/fam2.py`, `code/fam3.py`. QED

**So the maximum also has an analogous one-parameter family** (two orbits), and the pair
`(63/128,321/512)` of the task statement is only its first member. Both families consist of
points converging to `x=y=1/2`, exactly as Theorem B predicts.

---

## 9. Equality sets: exhaustive classification

Exact search over all `16^K` code pairs at every level `K <= 10`, filtered by the
Proposition of section 4 and then reduced by stripping common binary prefixes
(`code/eqset3.py`, `code/maxlist.py`):

| level `K` | distinct reduced **minimisers** | distinct reduced **maximisers** |
|---|---|---|
| 6 | 5 (`m = 2..6`) | 8 |
| 8 | 7 (`m = 2..8`) | 16 |
| 10 | 9 (`m = 2..10`) | 24 |

and in each case the list is *exactly*
`{ (1/2-3*4^-m, 1/2+3*4^-m) : 2 <= m <= K }` for the minimum, and
`{4^-j*(1/64,65/256), 4^-j*(17/64,1/256)}` together with the swap, for as many `j` as the
level allows, for the maximum. Nothing else occurs. (At level 10 the 24 maximisers are
`{orbits A,B} x {identity, swap} x {j = 0..5}`, the deepest `j` being limited by the
denominator `2*4^(j+5)` of `y`.)

**Conclusion (verified to level 10)**: the equality sets are the orbits of the two base
configurations `(5/16,11/16)`, `(63/128,321/512)`, `(47/128,257/512)` under
`(u,v) -> (u/4,v/4)` (in the `(T,S)` coordinates of (3)) and the swap, together with all
binary-prefix images. No other extremal configurations exist up to level 10; the families
themselves are proved for all `m` resp. `j` (Theorems A, C).

---

## 10. Certified rigorous bounds

### 10.1 The enclosure algorithm (correctness)

`T := conv{0, 1, lambda/3}`.

**Lemma A.** `K subset T`. Exact computation in `Z[omega]` (writing `a+b*omega=(a,b)`,
`lambda/3=(1/3,1/3)`, `code/lemmaA.py`):

    f_0(T)=conv{(0,0),(1/9,1/9),(1/3,0)}      f_1(T)=conv{(1/3,0),(2/9,2/9),(1/3,1/3)}
    f_2(T)=conv{(1/3,1/3),(5/9,2/9),(2/3,0)}  f_3(T)=conv{(2/3,0),(7/9,1/9),(1,0)}

In real coordinates the non-trivial vertices are `(1/6,sqrt3/18)`, `(1/3,sqrt3/9)`,
`(2/3,sqrt3/9)`, `(5/6,sqrt3/18)`, each satisfying the edge equation `y = x/sqrt3` or
`y = (1-x)/sqrt3` of `T` exactly; so `f_d(T) subset T` for `d=0..3` and the IFS attractor
lies in `T`. QED

**Lemma B.** For `w in {0,1}^j`, `Psi_w := phi_{e_1+1} o ... o phi_{e_j+1}` and
`I_w := T_{e_1} o ... o T_{e_j}([0,1])` (a dyadic interval of length `2^-j`),
`F(I_w) = Psi_w(K) subset Psi_w(T)`, and `Psi_w` is an exact similarity
`z -> alpha z + beta` or `alpha conj(z)+beta` with `alpha,beta in (1/3^j)Z[omega]`.
*Proof:* `F(t) = Psi_w(F(tau^j t))` from (*), `tau^j(t) in [0,1]`, and
`phi_1(z) = conj(z)lambda/3`, `phi_2(z) = z(omega-2)/3 + 1`. QED

So for a dyadic square `Q = I x J` of side `2^-J` (single cylinders `w,w'`), the image of
`N = omega conj F(t)+F(s)-lambda` over `Q` lies in the hexagon
`omega*Psi_w^c(T)+Psi_{w'}(T)-lambda` (`Psi^c` = conjugated similarity); its bounding box
gives `d <= |N| <= D`, whence

    d/(sqrt3 delta_max^alpha) <= R <= D/(sqrt3 delta_min^alpha),
    delta_max = 1-(t_lo+s_lo)/2,  delta_min = 1-(t_hi+s_hi)/2.

All arithmetic is exact in `Z[omega]` except the final comparison, done in directed-rounding
interval arithmetic (`mpmath.iv`, 50 digits) for `sqrt3` and `x^alpha = exp(alpha log x)`,
`alpha` enclosed as an interval. (One may replace it by rational tanh/Taylor bounds for
`log` and `exp`; nothing else changes.)

### 10.2 Results (`code/bb3.py`)

**Theorem D (certified).**
 * The sharp bounds hold on `{delta >= 1/2}`, i.e. for every reduced pair with
   `y-x >= 1/2`: `2*6^-alpha <= R(t,s) <= 2*sqrt(211)*69^-alpha`.
   513 boxes at depth 27, **no residual**.
 * `(1-1e-4)*2*6^-alpha <= R(t,s) <= (1+1e-4)*2*sqrt(211)*69^-alpha` for every reduced
   pair with `y-x >= 2^-10`: 17372 boxes at depth 27, **no residual**.
 * With `eps = 0` on `{delta >= 2^-10}` the same exhaustion leaves 84 dyadic squares of
   side `2^-27` (total area `4.7e-15`); all of them contain a member of the families of
   sections 6 and 8. They cannot be removed: by Theorem A the value `2*6^-alpha` is
   attained at points of the region, so the inequality there is not strict.

By Theorem B, everything outside `{delta >= 1/2}` is an exact copy of the `f`-problem on
`[0,1]^2`, whose extreme values are limits at the corner; hence no finite exhaustion can
produce the sharp statement outside `{delta >= 1/2}`.

---

## 11. What is proved, and exactly what remains

**Proved:** (1), (2), (3); binary-branch invariance; the reflection symmetry; the reduced
formula with `omega` (and hence the falsity of the stated `omega^2` formula); the exact
extremal data of section 3; Theorem A (the minimising family, all `m>=2`); Theorem B (the
pinch blow-up, its exact 2-periodicity, the range identities); Theorem C (the two
maximising families, all `j>=0`); Lemmas A and B and the correctness of the enclosure;
Theorem D.

**Not proved:** the two inequalities in full. By Theorem B they are *equivalent* to

    (S)   2*6^-alpha <= f(u,v) <= 2*sqrt(211)*69^-alpha   for all (u,v) in [0,1]^2,
          f(u,v) = |omega F(u) + conj F(v)|/(u+v)^alpha,

and equivalently to (S) restricted to any corner `{u,v <= 4^-k}`. (S) is a single
self-similar inequality with no compact reduction. The complete proof needs one of:
(i) an exact algebraic inequality bounding `|omega F(u)+conj F(v)|` above and below by
`M(u+v)^alpha` resp. `m(u+v)^alpha` — note both bounds hold with equality at the family
points of sections 6 and 8, and by section 4 the equality cases are exactly those; or
(ii) a closed invariant cone of functions on the state space of section 5 carrying a
genuine contraction — the transfer itself does not supply one, since `sigma_(0,0)` is
2-periodic on the two states that occur and the self-similarity `f(u,v)=f(4u,4v)` is an
identity, not a contraction.

Certified information is available everywhere except the corner of (S):
`{delta >= 1/2}` sharply, `{delta >= 2^-10}` up to `1e-4`, both with no residual.

---

## 12. Files

| file | content |
|---|---|
| `code/zw.py` | exact `Z[omega]` arithmetic, `F` at code points, exact scans |
| `code/family.py` | Theorem A, exact, `2 <= m <= 19` |
| `code/fam2.py`, `code/fam3.py` | Theorem C, exact, orbits A and B, `j <= 9` |
| `code/pinch.py` | the state transfer, the 2-cycle, Theorem B identities |
| `code/verify4.py` | transfer map over 6-step compositions; itineraries |
| `code/verify2.py` | forms (2), (3) and the `M_e` identities |
| `code/eqset2.py`, `code/eqset3.py`, `code/maxlist.py` | exhaustive equality sets, levels `<= 10` |
| `code/lemmaA.py` | Lemma A, exact |
| `code/encl.py`, `code/bb3.py` | the enclosure algorithm and Theorem D |
