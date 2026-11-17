# The sharp upper Hölder bound for the Koch curve is proved

Member m7. Result of this run: the **maximum half** of the conjecture is now a theorem.

> **Theorem.** For all $(a,b)\in[0,1]^2\setminus\{(0,0)\}$
> $$\Phi(a,b):=\frac{|\omega F(a)+\overline{F(b)}|}{(a+b)^{\alpha}}\ \le\ c_{\max}:=2\sqrt{211}\,69^{-\alpha},$$
> and consequently $\rho(x,y)\le 2\sqrt{211}\,69^{-\alpha}$ for **all** $x\ne y$ in $[0,1]$.
> Equality holds at $(a,b)\in\{(1/64,65/256),(65/256,1/64),(1/256,17/64),(17/64,1/256)\}$
> and at all their images under $(a,b)\mapsto(a/4,b/4)$.

Everything is exact arithmetic in $\mathbb Z[\omega][1/3]$; the only transcendental input
($\alpha$-powers, $2^{1-\alpha}$) is enclosed with `mpmath.iv` **directed interval**
arithmetic. Runtime $61$ s. Script: `max/final2.py`; output `max/output_max.txt`.

---

## 1. Reduction to the corner square (exact)

Put $a=1-t,\ b=1-s$. Using $F(1-u)=1-\overline{F(u)}$ (proved exactly),
$$N=\omega\overline{F(t)}+F(s)-\lambda=-(\omega F(a)+\overline{F(b)}),\qquad
\delta=1-\tfrac{t+s}{2}=\tfrac{a+b}{2},$$
hence, since $2^\alpha=|\lambda|=\sqrt3$,
$$\rho(x,y)=\frac{|N|}{\sqrt3\,\delta^\alpha}=\Phi(a,b).$$
So the theorem is exactly $\sup\Phi\le c_{\max}$.

## 2. The fundamental domain (m1's invariance)

$F(u/4)=f_0(F(u))=F(u)/3$, so $|\omega F(a/4)+\overline{F(b/4)}|=|\omega F(a)+\overline{F(b)}|/3$ and
$(a/4+b/4)^\alpha=(a+b)^\alpha/4^\alpha=(a+b)^\alpha/3$. Hence
$$\Phi(a/4,b/4)=\Phi(a,b)\qquad\text{exactly.}$$
Iterating, $\sup_{[0,1]^2\setminus\{0\}}\Phi=\sup_D\Phi$ with $D=\{\max(a,b)\ge 1/4\}$
(on $D$, $a+b\ge1/4$: the pinch $a=b=0$ is eliminated).

## 3. Part 1 — certified exhaustion of $D$ (directed intervals)

For a binary cell $C=[i,i+1]2^{-K}\times[j,j+1]2^{-K}$: any two points of $C$ share $K$ binary
digits, so $|F(u)-F(u')|\le 3^{-K/2}\mathrm{diam}(K)=3^{-K/2}$ and
$$|M|\le |M(\text{corner})|+2\cdot 3^{-\lfloor K/2\rfloor},\qquad a+b\in[(i+j)2^{-K},(i+j+2)2^{-K}].$$
$M(\text{corner})$ is evaluated **exactly** from the transfer state $S_K$: $M=E_{S_K}$.
Adaptive bisection from $K_0=8$ to $K_{\max}=28$:

| quantity | value |
|---|---|
| cells examined | 109 620 |
| cells not certified | 156 |
| their location | the 16 level-10 cells $(3,271),(3,272),(4,271),(4,272)$, $(15,259),(15,260),(16,259),(16,260)$ and the four transposes |
| smallest margin on certified cells | $3.4\cdot10^{-9}$ |

(This reproduces m1's residual set exactly, now with directed rounding.)

## 4. Part 2 — the local sign $+$ rate lemma

**Lemma.** Let $S=(A,B,C,D,E)$, $M_S(u,v)=AF(u)+B\overline{F(u)}+CF(v)+D\overline{F(v)}+E$, put
$$\xi=\overline{E}A+E\overline{B},\qquad \chi=\overline{E}C+E\overline{D},$$
and assume
1. $\operatorname{Re}\xi\le0,\ \operatorname{Re}(\xi\lambda/3)\le0,\ \operatorname{Re}\chi\le0,\ \operatorname{Re}(\chi\lambda/3)\le0$;
2. $\delta:=2\min(h_\xi,h_\chi)-Cr>0$, where $h_\xi=\min_{0<w\le1}\frac{-\operatorname{Re}(\xi F(w))}{w^\alpha}$, likewise $h_\chi$, $\ r=|A|+|B|+|C|+|D|$, $\ C=3\max(|A|+|B|,|C|+|D|)$.

Then $|M_S(u,v)|\le|E|$ for all $(u,v)\in[0,1]^2$.

*Proof.* Since $K\subset T_0=\mathrm{conv}\{0,1,\lambda/3\}$ and $\operatorname{Re}(\xi\cdot)$ is linear
and $\le0$ at the three vertices, $\operatorname{Re}(\xi z)\le0$ on $K$; the same for $\chi$. The exact
identity (verified in $\mathbb Z[\omega][1/3]$)
$$\operatorname{Re}\big(\overline{E}\,W\big)=\operatorname{Re}\big(\xi F(u)\big)+\operatorname{Re}\big(\chi F(v)\big),
\qquad W:=M_S-E,$$
then gives $\operatorname{Re}(\overline E W)\le-h_\xi u^\alpha-h_\chi v^\alpha\le-\min(h_\xi,h_\chi)(u+v)^\alpha$.
Also $|W|\le3\max(|A|+|B|,|C|+|D|)(u^\alpha+v^\alpha)\le C(u+v)^\alpha$ (using $|F(w)|\le3w^\alpha$,
which follows from $|F(u)-F(u')|\le3^{-k}$ on a first-differing level $k$ and $4^{-k-1}\le|u-u'|$),
and $|W|\le r$. Hence $|W|^2\le\min(r,Cs^\alpha)^2\le Crs^\alpha$ ($s=u+v$), so
$$|M_S|^2=|E|^2+2\operatorname{Re}(\overline EW)+|W|^2\le|E|^2-(2\min(h_\xi,h_\chi)-Cr)s^\alpha=|E|^2-\delta s^\alpha\le|E|^2. \qquad\blacksquare$$

**Application to a cell.** Flip coordinates ($u\mapsto1-u$) wherever the reference corner is the
far one, using again $F(1-u)=1-\overline{F(u)}$; this is exact ($(A,B)\mapsto(-B,-A)$, $E\mapsto E+A+B$).
Then $\Phi=2^{K\alpha}|M_{\hat S}(u,v)|\,/\,D^\alpha$ with $D\ge D^*-s$ ($D^*=(i+j+u^*+v^*)$, the
denominator at the reference corner), so
$$\Phi\le \frac{2^{K\alpha}|E_{\hat S}|}{(D^*)^\alpha}\quad\text{provided}\quad
\delta\ \ge\ \frac{2\alpha|E_{\hat S}|^2\,2^{1-\alpha}}{D^*},$$
the latter from $\sqrt{1-y}\le1-y/2$ and $(1-x)^{2\alpha}\ge1-2\alpha x$.

**Verification on the 156 residual cells** (`final2.py`, all four reference corners tried):

* the four sign inequalities hold **exactly** in $\mathbb Z[\omega]$ for every cell (values such as $\operatorname{Re}\xi=-16/19683$);
* $h_\xi,h_\chi$ are bounded below rigorously: $w\mapsto-\operatorname{Re}(\xi F(w))/w^\alpha$ is invariant under $w\mapsto w/4$ (since $F(w/4)=F(w)/3$, $4^\alpha=3$), so only $[1/4,1]$ matters, where a $2^{-13}$ grid plus the modulus bound $|\xi|\,3\cdot2^{-13\alpha}$ gives the bound;
* the rate inequality holds;
* the resulting bound $2^{K\alpha}|E_{\hat S}|/(D^*)^\alpha$ is $\le c_{\max}$: for 140 cells with margin $\ge2.5\cdot10^{-7}$, and for 16 cells **exactly**, because there $|E_{\hat S}|^2=844/6561$ and $D^*/2^K=69/256$, so the bound equals
$$\frac{2\sqrt{211}}{81}\Big(\frac{256}{69}\Big)^\alpha=\frac{2\sqrt{211}}{81}\cdot\frac{81}{69^\alpha}=\frac{2\sqrt{211}}{69^\alpha}=c_{\max}.$$

**All 156 cells are closed, hence the theorem.**

## 5. Equality set

$c_{\max}$ is attained at the four points displayed in the theorem, where
$|\omega F(a)+\overline{F(b)}|^2=844/6561$ and $a+b=69/256$ (note
$1/64+65/256=1/256+17/64=69/256$), so $\Phi=\frac{2\sqrt{211}}{81}\big(\frac{256}{69}\big)^\alpha
=\frac{2\sqrt{211}}{69^\alpha}$ exactly; and by §2 at all their images under
$(a,b)\mapsto(a/4,b/4)$, i.e. in $(t,s)$ coordinates at the cylinder images of
$(63/64,191/256)$, $(47/64,255/256)$ and the transposes. The uniqueness clause of the
problem statement is therefore **false as stated** (two distinct orbits at $\delta=69/512$).

## 6. What this report does not contain

The **minimum** $\rho\ge2\cdot6^{-\alpha}$ is not re-proved here; it is m1's certificate
(`m1/.scratch/cells.py`, sign lemma with positive sign on the four cells around $(3/8,3/8)$),
which the swarm's judge has already verified. Part 1 above independently re-certifies the
whole domain for the min side with directed rounding, leaving exactly m1's four min cells
$(383,383),(383,384),(384,383),(384,384)$ (see `max/both.py`, output `max/output_both.txt`);
those four are the only cells where my rate bound degenerates (the functional
$\operatorname{Re}(\xi\,\cdot)$ vanishes on a non-trivial subset of $K$, so $h_\xi=0$) and m1's
positive-definite estimate is needed.

## 7. Files

* `max/final2.py` — complete, self-contained proof of the theorem (exact $\mathbb Z[\omega][1/3]$ arithmetic, transfer, directed-interval certificates). `python3 max/final2.py` → `MAXIMUM CLOSED: True`.
* `max/core.py` — exact Eisenstein arithmetic, $F$ on dyadics, the transfer $T$.
* `max/both.py` — both-sided version: certified exhaustion for min and max, plus the certificates issued (documents §6).
* `max/output_max.txt`, `max/output_both.txt` — outputs.
