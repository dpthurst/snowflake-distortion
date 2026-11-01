# Related prior work for *Expansion*

*Prepared by Wren, an Open Math OS Agent.*

**Scope.** I read the draft at `dpthurst/expansion`, commit `73fe2f1` (2026-10-18), including `alternatives.tex` and the unintegrated maximum-growth material in `clippings.tex`. This memo is a literature map for the two threads you named. Every work it cites is collected once, by thread, in the bibliography at the end (Section 9). One distinction runs through the whole report and is easy to lose: **topological** invariants (invariant under topological conjugacy — what your paper produces), **metric-bound** invariants (invariant only under bi-Lipschitz changes of the metric), and **smooth** invariants (which need a differentiable structure). Section 2.1 sets out the three kinds and places the literature's quantities in them; each later subsection says which kind it is discussing. It deliberately does not address any open question in the draft, and it does not try to prove or disprove anything in it.

**One-paragraph summary.** Both threads already have substantial literatures, and the two are more connected than the draft's separation suggests. As far as I found, the closest published invariant to your $h_{\min}$ is the exponential decay rate of **Lebesgue numbers** of the refined covers (Sun 2012): your in-radius $r(\mathcal U;\mathcal U^{(n)})$ is essentially the reciprocal of the Lebesgue number $\delta(\mathcal U^{(n)})$, so, up to the comparison between in-radius and Lebesgue number, $h_{\min}$ is minus Sun's $h_L$. The closest smooth invariant to your $\lambda_{\exp}$ is **Hirsch's logarithmic expansion rate** $E(f,K)=\sup_n \frac1n\log EC(f^n,K)$ and **Schreiber's theorem** that, for $C^1$ maps, $E$ is the infimum of the smallest Lyapunov exponent over ergodic invariant measures. The closest circle of ideas to your visual-metric theorem is the **entropy–dimension–expansion inequality** (Fathi 1989; Dai–Zhou–Geng 1998; Misiurewicz 2004; Fujita–Kato–Matsumoto 2010; Roth–Roth 2020). In its sharp form it says that for a metric in which $f$ expands small distances by exactly $\lambda$,
$$
\dim X\cdot\log\lambda \;\le\; D_d(X)\log\lambda \;\le\; h_{\mathrm{top}}(f) \;\le\; D_d(X)\log\operatorname{Lip}_d(f),
$$
where $D_d(X)$ is the box dimension in the metric $d$. Your $h_{\min}$ is the topological minimal distance growth, $h_{\mathrm{top}}$ is the volume growth, and $\lambda_{\exp}=e^{h_{\min}}$ is the maximal expanding constant. The draft's theorem $\log\lambda_{\exp}=h_{\min}$ is the topological counterpart of Schreiber's smooth theorem. At the other end, the closest published objects are Kifer's pointwise maximum and minimum characteristic exponents in metric spaces (1983) — the pointwise analogue of your $h_{\max}$, but depending on the metric — and the infimum of the Lipschitz constant in a conjugacy class of interval maps (Bobok–Roth 2019), the closest relative of $\lambda_{\mathrm{Lip}}$; Section 8 maps that side, and nothing I found asserts or denies $h_{\max}=\log\lambda_{\mathrm{Lip}}$.

---

## 1. The six papers I would read first

1. **Peng Sun, "Exponential decay of Lebesgue numbers,"** *Discrete Contin. Dyn. Syst.* **32** (2012), 3773–3785; arXiv:1002.0549. Defines
   $$
   h_L^\pm(f,\mathcal U)=\limsup/\liminf_{n\to\infty}-\frac1n\log\delta\bigl(\textstyle\bigvee_{k<n}f^{-k}\mathcal U\bigr),
   $$
   where $\delta$ is the Lebesgue number, and $h_L^\pm(f)=\sup_{\mathcal U}h_L^\pm(f,\mathcal U)$. Proves $\dim_B(X)\,h_L^-(f)\ge h(f)$ and $\dim_H(X)\,h_L^+(f)\ge h(f)$; for Lipschitz $f$, $h_L^+\le\log L(f)$; $h_L$ is bi-Lipschitz invariant but not a topological invariant. This is the closest published invariant to your $h_{\min}$.

2. **Sebastian J. Schreiber, "Expansion rates and Lyapunov exponents,"** *Discrete Contin. Dyn. Syst.* **3** (1997), 433–438; DOI 10.3934/dcds.1997.3.433. For a $C^1$ endomorphism $f$ of a Riemannian manifold and a compact positively invariant $K$ on which $f$ is immersive, defines the expansion constant $EC(f,K)$ (the largest $\rho\ge0$ such that for every $0<\rho^*<\rho$ there is $\delta>0$ with $d(fx,fy)\ge\rho^*d(x,y)$ whenever $x\in K$, $y\in M$ and $d(x,y)\le\delta$) and Hirsch's logarithmic expansion rate $E(f,K)=\sup_n\frac1n\log EC(f^n,K)$. Proves
   $$
   E(f,K)=\inf_\mu E(f,\mu),
   $$
   the infimum of the smallest Lyapunov exponent over ergodic invariant measures supported on $K$. This is the smooth analogue of $\log\lambda_{\exp}=h_{\min}$.

3. **Chikara Fujita, Hisao Kato, Masahiro Matsumoto, "Fractal metrics of Ruelle expanding maps and expanding ratios,"** *Topology Appl.* **157** (2010), 615–628; DOI 10.1016/j.topol.2009.11.002. For a Ruelle expanding map $f$ of a compactum $X$ and every $s>1$, constructs an admissible metric $d$ and an **exact** expansion ratio $\lambda\in(1,s)$: $d(fx,fy)=\lambda d(x,y)$ for $d(x,y)\le\epsilon$. In that metric,
   $$
   \dim_H(X,d)=\underline D_d(X)=D_d(X)=\frac{h(f)}{\log\lambda}.
   $$
   Proposition 4.1 gives, for $1<\lambda_2\le\lambda_1$ and $\lambda_2 d(x,y)\le d(fx,fy)\le\lambda_1 d(x,y)$ for small $d(x,y)$,
   $$
   D_d(X)\log\lambda_2\le h(f)\le D_d(X)\log\lambda_1.
   $$
   Remark 4.4 deduces $\log\lambda\le h(f)/\dim X$ for a strictly expanding ratio $\lambda$; hence the set of expanding ratios is bounded above by $e^{h/\dim X}$. This is the closest "set of possible expanding constants" result.

4. **Albert Fathi, "Expansiveness, hyperbolicity and Hausdorff dimension,"** *Comm. Math. Phys.* **126** (1989), 249–262; DOI 10.1007/BF02125125. For a compact hyperbolic set $K$ of a $C^1$ diffeomorphism, with uniform contraction rate $\rho<1$ for $f$ on the stable manifolds and for $f^{-1}$ on the unstable manifolds,
   $$
   \mathrm{HD}(K)\le \frac{2\,\mathrm{ent}(f|_K)}{\log(1/\rho)}.
   $$
   The same argument with no stable direction gives, for an expanding map, $\mathrm{HD}(K)\le h(f)/\log\lambda$. The paper also constructs a "hyperbolic metric" and uses the inequality to prove Mañé's theorem that expansive homeomorphisms have finite topological dimension.

5. **Samuel Roth, Zuzana Roth, "Inequalities for entropy, Hausdorff dimension, and Lipschitz constants,"** *Studia Math.* **250** (2020), 253–264; arXiv:1807.01196. Defines the HausLip constant
   $$
   \mathrm{HausLip}(X,f)=\inf_{d\in\mathcal D(X)}\mathrm{HD}_d(X)\log^+\operatorname{Lip}_d(f),
   $$
   proves it is a topological conjugacy invariant and $\ge h(f)$, and proves $\mathrm{HausLip}(X,f)=h(f)$ for positively expansive maps and for linear maps of tori. The proof for expansive maps constructs a metric in which $\operatorname{Lip}_D(f^n)\le16\,\mathrm{SkewLoc}_D(f^n)$, i.e. a metric that is uniformly expanding up to a bounded factor. This is the closest published construction to your visual-metric theorem.

6. **S. A. Dovbysh, "Optimal Lyapunov metrics of expansive homeomorphisms,"** *Izv. Math.* **70** (2006), 883–929; DOI 10.1070/IM2006v070n05ABEH002332. Sharpens Reddy, Sakai, and Fried: every expansive homeomorphism of a metrizable compactum admits a compatible Lyapunov metric; for A$^\#$-homeomorphisms the lower bounds for the contraction and expansion constants are attained simultaneously by one metric satisfying Ruelle's Lipschitz axiom. This is the two-sided analogue of optimizing the metric in your $\lambda_{\exp}$.

---

## 2. Your objects in the literature's language

| Draft object | Closest literature object | Main references |
|---|---|---|
| $(\lambda,\epsilon)$-expanding | expanding in the sense of Ruelle; distance-expanding; positively expansive + expanding metric | Reddy 1982; Przytycki–Urbański 2010; Shub 1969 |
| $\lambda_{\exp}(f)$ | maximal expansion constant; logarithmic expansion rate $E$; expanding ratio | Hirsch 1994; Schreiber 1997; Fujita–Kato–Matsumoto 2010; Dovbysh 2006 |
| $d(x,y;\mathcal U)$ | chain distance in a cover; uniformity chain | Bowen 1971/1973; Misiurewicz 1986; uniform-space entropy |
| $h_{\mathrm{dist}}(x,y)$ | pairwise Lyapunov exponent / pointwise expansion rate | Oseledets; Hirsch 1994; Schreiber 1997 |
| $r(\mathcal U;\mathcal V)$ | reciprocal of the Lebesgue number | Sun 2012 |
| $h_{\min}$ | minus the Lebesgue-number decay rate; infimum Lyapunov exponent (smooth case) | Sun 2012; Schreiber 1997 |
| $h_{\mathrm{diam}}$, $h_{\max}$ | diameter growth; maximal characteristic exponent; $R(f)=\lim\frac1n\log\sup_x\lVert Df^n_x\rVert$ | Kifer 1983; Newhouse 1988; Yomdin 1987; Pacifico–Vieitez 2020 |
| $\lambda_{\mathrm{Lip}}$ | infimum of Lipschitz constants over length metrics; conjugacy-class Lipschitz infimum; HausLip constant | Bobok–Roth 2019; Roth–Roth 2020; Milnor §13 |
| $h_{\mathrm{top}}$ | topological entropy; volume growth | Adler–Konheim–McAndrew 1965; Bowen 1971; Newhouse 1988; Yomdin 1987 |

One more caution on the table: the draft's $\lambda_{\mathrm{Lip}}$ is an infimum over *length* metrics, while the Bobok–Roth and Roth–Roth literature usually optimizes over all compatible metrics or over Euclidean metrics in a conjugacy class. The quantities are closely related but not literally the same; the table is a map of ideas, not an identity of definitions. Section 8 takes the maximum-end rows up in detail, including why the length-metric restriction is not optional.

Three terminology warnings.

- **"Expansive" vs "positively expansive."** In the literature, *expansive* usually means the two-sided homeomorphism notion: $d(f^n x,f^n y)<\gamma$ for all $n\in\mathbb Z$ implies $x=y$. Your $\mathbb N$-expansive is *positively expansive*. For a homeomorphism of a compact space, positive expansiveness is extremely restrictive (it forces a finite space; see the note "Positively expansive homeomorphisms of compact spaces," arXiv:math/0308294, and the references there). Your one-sided setting is the right one for non-invertible maps, but the word "expansive" alone will be read as the two-sided notion by many dynamicists; the draft already writes "$\mathbb N$-expansive or just expansive," which is a good guard.
- **"Expansive constant" vs "expansion constant."** An *expansive constant* is usually the **scale** $\epsilon$ (Sun 2011, Walters, Gottschalk–Hedlund). An *expansion constant* / *expanding ratio* / *expansion rate* is the **rate** $\lambda$ (Hirsch, Schreiber, Fujita–Kato–Matsumoto). Your $\lambda_{\exp}$ is a rate; your $\epsilon$ is a scale. The phrase "maximal expanding constant" is ambiguous between the two, and the literature uses both.

- **The refinement order runs the other way in your draft.** You write $\mathcal U\le\mathcal V$ when $\mathcal U$ refines $\mathcal V$, so your $\vee$ is the common coarsening and your $\wedge$ the common refinement. In dynamics the opposite order is universal: refinement makes a cover *larger*, the common refinement is the join, and $N$, $H$ and $h(\cdot;\mathcal U)$ all increase under refinement. Adler–Konheim–McAndrew write $\mathfrak{A}\prec\mathfrak{B}$ for "$\mathfrak{B}$ refines $\mathfrak{A}$", and their Property 3 fixes the direction without any guesswork: they prove $\mathfrak{A}<\mathfrak{B}\Rightarrow H(\mathfrak{A}\vee\mathfrak{B})=H(\mathfrak{B})$, and the common refinement has the entropy of the *finer* cover, so $\mathfrak{B}$ is the finer one and $<$ points from coarse to fine. Sun (2012) writes $\mathcal U\succ\mathcal V$ for "$\mathcal U$ refines $\mathcal V$" (his Lemma 2.2: refinement decreases the Lebesgue number); Bruin's notes write $\mathcal U\succeq\mathcal V$; Walkden's notes write $\alpha\le\beta$ for "$\beta$ refines $\alpha$"; the mean-dimension literature writes $\beta\succ\alpha$ for "$\beta$ refines $\alpha$". The order with refinement *decreasing* a cover, and the common refinement as the meet, is the one used in covering uniformities (Willard; Isbell: a uniformity is closed under coarsening, and its meet is the common refinement) — not the one used in dynamics. The notion of refinement itself is of course the same everywhere (Engelking, §1.5; Willard, §36); only its order differs. Two practical consequences. First, the glyph does not identify the order: $\prec$ means "refines" in general topology and "is coarser than" in Adler–Konheim–McAndrew, and $\le$ appears in both directions in the sources above, so a paper has to say which order it uses. Second, once the order is declared the rest follows from the standard lattice meanings — $\wedge$ is the meet, $\vee$ the join — so nothing further needs saying, and your draft does declare it: `covers.tex` writes $\mathcal V\le\mathcal U$ for "$\mathcal V$ refines $\mathcal U$", and its remark after the definition of $\mathcal U\wedge\mathcal V$ records the reversal from AKM together with the corresponding reversal of meet and join. The reader's whole task is then to apply the lattice meanings to that order: the common refinement, which dynamics writes $\vee$, is the *meet* here and appears as $\wedge$ — it is the operation in $\mathcal U_n=\mathcal U\wedge f^{-1}\mathcal U\wedge\cdots\wedge f^{-(n-1)}\mathcal U$, and $\vee$ itself never occurs in the paper outside that remark. The convention is also consistent at the empty end: $\mathcal U_0=\{X\}$ is the coarsest cover, which is what an empty meet should be in this order.

**Where the two orders come from.** Adler–Konheim–McAndrew's footnote to "reflexive partial ordering" reads "For terminology see J. L. Kelley, *General topology*, Van Nostrand, New York, 1955", so they present the terminology as inherited. But Kelley attaches no inequality to refinement of *covers* — he states the relation in words — and his own covering-uniformity condition (Problem J, p. 207) is closed under coarsening, which is the orientation topology kept. The only refinement order Kelley does symbolize is for subdivisions in the chapter on integration, where $S\ge S'$ iff $S$ refines $S'$: the direction AKM use. What AKM's direction matches is the *partition* convention of ergodic theory, and their paper says so: the common refinement of two covers is called the join (Definition 2), the remark after Property 1 has $\mathfrak{A}<\mathfrak{A}\vee\mathfrak{B}$ and $\mathfrak{B}<\mathfrak{A}\vee\mathfrak{B}$ and concludes that the covers form a directed set, and in §4 the same $\vee$ is used for measurable partitions, "the common refinement", with $H_\mu(\mathfrak{A}\vee\mathfrak{B})\le H_\mu(\mathfrak{A})+H_\mu(\mathfrak{B})$. So the dynamics order is AKM's own alignment of covers with partitions; it is not a misreading of Kelley, and it is not the covering-uniformity order. The partition convention itself is older than the covers paper — it comes from the founding papers of metric entropy, Kolmogorov's 1958 note and Sinai's 1959 note — and is written out explicitly in Rokhlin's survey (§1.3): he fixes $\xi<\xi'$ for "$\xi'$ is a subpartition of $\xi$", defines the product $\textstyle\bigvee_\alpha\xi_\alpha$ as the join (every $\xi_\alpha$ is below it, and it is the least such partition) and the intersection $\textstyle\bigwedge_\alpha\xi_\alpha$ as the meet, with the point partition $\varepsilon$ at the top and the trivial partition $\nu$ at the bottom. So the measure-theoretic convention and AKM's cover convention were already the same one when AKM wrote: their paper was received in August 1963, and its first definition says it is conforming to ergodic theory. (The two Doklady notes are not on mathnet and are behind a paywall, so I have not read them; the earliest statement in print I have seen is Rokhlin's.)

### 2.1 Three kinds of invariant, and which is which

Three kinds of object appear in this report, and they should not be conflated:

| Kind | Invariance | Examples in this report |
|---|---|---|
| **Topological** | invariant under topological conjugacy | $h_{\mathrm{top}}$; your $h_{\min}$, $h_{\mathrm{dist}}$, $h_{\mathrm{diam}}$, $h_{\max}$; your $\lambda_{\exp}$ and $\lambda_{\mathrm{Lip}}$; $\mathrm{HausLip}$ (Roth–Roth); $\Lambda(f)$ (Bobok–Roth); $\lambda(f)$ (Braga–Sieg) |
| **Metric-bound** | invariant under bi-Lipschitz changes of the metric, not under topological conjugacy | Sun's $h_L^\pm$ and $h_E^\pm$; Kifer's $\Lambda^\pm(x)$, $\lambda^\pm(x)$ in a fixed metric; Morales–Thieullen–Villavicencio's $\chi_f$; Hirsch's one-step $EC(f,K)$ |
| **Smooth** | needs a differentiable structure; independent of the choice of Riemannian metric | the Lyapunov exponents $\lambda_i$; Schreiber's $E(f,K)$; Newhouse–Yomdin's $v_k(f)$ and $R(f)$; the Lyapunov spectrum $[\alpha^-,\alpha^+]$; Bowen's equation and the pressure formula |

**Every invariant in your draft is topological**, and the reason is the same in each case: it is built either from open covers (the $h$'s and $h_{\mathrm{dist}}$) or as a supremum or infimum over *all* compatible metrics ($\lambda_{\exp}$, $\lambda_{\mathrm{Lip}}$). A homeomorphism carries the class of compatible metrics to itself, so a supremum over that class is automatically invariant under topological conjugacy. That is a design decision, not an accident, and it is what separates $\lambda_{\exp}$ from every smooth expansion rate in Section 4. The distinction is precisely *all compatible metrics* versus *all Riemannian metrics*: the first supremum is topological, the second is a smooth invariant — Schreiber's $E$ (Section 4.4). The two are suprema over different classes, but the inclusion does give one inequality, as you pointed out: a Riemannian metric in which $f^n$ expands by $\rho$ yields a *compatible* metric in which $f$ expands by $\rho^{1/n}$, so $\lambda_{\exp}\ge e^{E}$ (the one-line proof is in Section 4.4). The reverse inequality is the open half; by your theorem it is the same as $h_{\min}\le E$.

**Most of the literature's quantities are not topological.** Sun's $h_L$ is invariant under bi-Lipschitz changes of the metric and not under topological conjugacy (his Examples 6.4–6.5); his $h_E$ is metric-bound too, since it is read off a separation scale in a fixed metric; Kifer's exponents and Morales–Thieullen–Villavicencio's $\chi_f$ are bi-Lipschitz invariant only (Section 8.4); the Lyapunov exponents, Schreiber's $E$, and Newhouse–Yomdin's $v_k$ and $R$ need a differentiable structure. So "your $h_{\min}$ is Sun's $h_L$" is a comparison between a topological invariant and a metric-bound one, and the honest statement is the comparison $h_{\min}\approx-h_L$ *in a fixed metric and for an expansive cover*, not an identity of invariants.

**The interesting theorems are the bridges between the columns.**

- *smooth to smooth*: Schreiber's $E(f,K)=\inf_\mu E(f,\mu)$ identifies Hirsch's asymptotic metric expansion rate with the infimum of the Lyapunov exponents.
- *topological to metric-bound*: Reddy's positively expansive $\Leftrightarrow$ expands small distances in some metric, and Fathi's hyperbolic metric for expansive homeomorphisms.
- *topological to metric*: the entropy–dimension–expansion inequalities (Fathi, Dai–Zhou–Geng, Misiurewicz, Fujita–Kato–Matsumoto, Roth–Roth), which relate the topological $h_{\mathrm{top}}$ to a dimension in a chosen metric times a logarithmic expansion or Lipschitz rate.
- *topological to topological*: your theorem $\log\lambda_{\exp}=h_{\min}$ — both sides are topological, and the content is that a supremum over metrics equals a cover-chain growth rate.

Each subsection below says which kind it is discussing.

---

## 3. Theme 1: entropy-like invariants built from distances

### 3.1 Entropy from separated and spanning sets (Bowen, Dinaburg)

The standard distance-based definition of topological entropy is
$$
d_n^f(x,y)=\max_{0\le k<n}d(f^k x,f^k y),\qquad
h_{\mathrm{top}}(f)=\lim_{\epsilon\to0}\limsup_{n\to\infty}\frac1n\log N(n,\epsilon),
$$
where $N(n,\epsilon)$ is the maximal cardinality of an $(n,\epsilon)$-separated set (or the minimal cardinality of an $(n,\epsilon)$-spanning set). This is "distance" in the sense that it uses the metric, but it is still a **volume** count: it counts how many $\epsilon$-separated orbit segments there are. The metric enters the definition, but the $\epsilon\to0$ limit removes it: $h_{\mathrm{top}}$ is topological, and this is the one place in Theme 1 where a metric-based definition and a topological outcome agree exactly.

Your $d(x,y;\mathcal U^{(n)})$ is the cover analogue of $d_n^f(x,y)$, and your $h_{\mathrm{dist}}(x,y;\mathcal U)$ is the exponential growth rate of a **single pair's** $d_n^f$-distance, not the count of separated sets. The standard entropy is the sup over scales of the growth rate of the number of pairs; your $h_{\mathrm{dist}}$ is the growth rate of one pair. The two are related by the variational principle and by the Lyapunov spectrum (see §3.7), but they are not the same construction.

### 3.2 Entropy from $\epsilon$-chains (Misiurewicz 1986)

A genuinely chain-based definition of entropy is due to Misiurewicz. For $\delta,\epsilon>0$, let $S(n,\delta,\epsilon)$ be the maximal cardinality of a set of $\delta$-separated $\epsilon$-chains of length $n$ (an $\epsilon$-chain is a sequence $x_0,\dots,x_{n-1}$ with $d(fx_i,x_{i+1})\le\epsilon$; two chains are $\delta$-separated if $d(x_i,y_i)\ge\delta$ for some $i$). Then
$$
h(f)=\lim_{\delta\to0}\lim_{\epsilon\to0}\limsup_{n\to\infty}\frac1n\log S(n,\delta,\epsilon).
$$
This is the closest classical "chain" definition to your cover-chain construction, but it counts chains rather than tracking the chain distance of one pair. It shows that chain distances can recover the volume entropy, which is a useful sanity check for any chain-based definition. A modern account with this formula is Richeson–Wiseman (2008). The formula is topological as well: the $\delta\to0$, $\epsilon\to0$ limits remove the metric, so this is a chain formula for the topological invariant, not a metric-bound quantity.

### 3.3 Entropy from Lebesgue numbers (Sun 2012): the closest invariant to $h_{\min}$

For an open cover $\mathcal U$ of a compact metric space, let $\delta(\mathcal U)$ be its Lebesgue number: the largest $\delta>0$ such that every ball of radius $\le\delta$ is contained in some element of $\mathcal U$. Put $\delta_n(f,\mathcal U)=\delta(\bigvee_{k<n}f^{-k}\mathcal U)$, and define
$$
h_L^-(f,\mathcal U)=\liminf_{n\to\infty}-\frac1n\log\delta_n(f,\mathcal U),\qquad
h_L^+(f,\mathcal U)=\limsup_{n\to\infty}-\frac1n\log\delta_n(f,\mathcal U),
$$
and $h_L^\pm(f)=\sup_{\mathcal U}h_L^\pm(f,\mathcal U)$ over finite open covers.

Sun proves:

- $h_L^\pm(f)\ge0$, $h_L^-(f)\le h_L^+(f)$;
- $h_L^\pm(f^n)=n\,h_L^\pm(f)$;
- if $\mathcal U$ refines $\mathcal V$ (in Sun's notation $\mathcal U\succ\mathcal V$), then $h_L^\pm(f,\mathcal U)\ge h_L^\pm(f,\mathcal V)$;
- $\dim_B(X)\,h_L^-(f)\ge h(f)$ and $\dim_H(X)\,h_L^+(f)\ge h(f)$;
- if $f$ is Lipschitz with constant $L(f)$, then $h_L^+(f,\mathcal U)\le\log L(f)$ for every $\mathcal U$;
- $h_L^\pm$ is invariant under bi-Lipschitz conjugacy, but **not** a topological invariant (Examples 6.4–6.5 give topologically conjugate homeomorphisms with different $h_L$);
- if $f$ is an expansive homeomorphism with expansive constant $\gamma$, then every open cover of diameter $<\gamma$ computes $h_L^\pm(f)$ (Corollary 3.9).

**Comparison with your in-radius.** Your $r(\mathcal U;\mathcal V)$ is the largest $k$ such that every $\mathcal V$-chain of length $\le k$ lies in a single $\mathcal U$-element. If $\mathcal V$ has mesh $\eta$ and $\mathcal U$ has Lebesgue number $\delta$, then a $\mathcal V$-chain of length $k$ has diameter at most $k\eta$, so $k\eta<\delta$ implies $r(\mathcal U;\mathcal V)\ge k$; conversely, a maximal failing chain gives a lower bound of the same shape. Thus
$$
r(\mathcal U;\mathcal V)\asymp \frac{\delta(\mathcal U)}{\operatorname{mesh}(\mathcal V)}
$$
up to constants, and for the dynamical sequence $\mathcal V=\mathcal U^{(n)}$,
$$
h_{\min}=\lim_{n\to\infty}\frac1n\log r_{\mathcal U}(n)
\;\approx\;-\lim_{n\to\infty}\frac1n\log\delta_n(f,\mathcal U)
=-h_L(f,\mathcal U).
$$
The approximation is not an identity: it uses that, for an expansive cover in a compatible metric, the mesh and the Lebesgue number of $\mathcal U^{(n)}$ decay at the same exponential rate. That is exactly the situation in which your theorem makes $h_{\min}$ independent of $\mathcal U$. Sun's Corollary 3.9 is the analogous statement for $h_L$ and expansive homeomorphisms.

This is, as far as I found, the closest prior work to your $h_{\min}$. The main difference is that Sun works with a fixed metric and a one-cover Lebesgue number, and his $h_L$ is **not** a topological invariant but a metric-bound one, invariant under bi-Lipschitz changes of the metric only; your $h_{\min}$ is defined by minimizing over pairs and suping over covers, and your theorem says it *is* a topological invariant for expansive maps. In the language of the draft, Sun's $h_L$ is the metric-bound upper/lower version of $-\log\lambda_{\exp}$; passing from $h_L$ to $h_{\min}$ is exactly the passage from a metric invariant to a topological one.

### 3.4 Entropy from expansive constants (Sun 2011)

A related but distinct "distance" invariant is the exponential decay rate of the **expansive constant**. For an expansive map $f$, let $\gamma(f)$ be the largest $\epsilon>0$ such that $d(f^n x,f^n y)<\epsilon$ for all $n\ge0$ implies $x=y$. Define
$$
h_E^+(f)=\limsup_{n\to\infty}-\frac1n\log\gamma(f^n),\qquad
h_E^-(f)=\liminf_{n\to\infty}-\frac1n\log\gamma(f^n).
$$
Sun proves $h_E^\pm(f)\ge0$, $h_E^-(f)\le h_E^+(f)$, $h_E^\pm(f)=h_E^\pm(f^{-1})$ for homeomorphisms, $h_E^+(f^n)\le n h_E^+(f)$, $h_E^-(f^n)\ge n h_E^-(f)$, $h_E^\pm(f)\le\log L(f)$ for Lipschitz $f$, and the dimension bound $h_E^-(f)\,\dim_B^+X\ge h(f)$, $h_E^+(f)\,\dim_B^-X\ge h(f)$. For a two-sided subshift of finite type with metric $d(\omega,\omega')=q^{-\min\{|j|:\omega_j\ne\omega'_j\}}$, $h_E^\pm(\sigma_A)=\frac12\log q$, so $h_E\cdot\dim_H=h_{\mathrm{top}}$; the one-sided case is analogous with the corresponding dimension.

The expansive constant is a **separation scale**, not an in-radius: it is the largest scale at which distinct orbits are eventually separated, whereas your $r_{\mathcal U}(n)$ is the largest number of steps that fit inside a cover element. The two are related (both measure how fast the dynamics resolves the space), but they are not the same invariant. Sun's inequality $h_E\le h_L$ is a relation between the two distance-based quantities. Both are metric-bound, not topological: each is read off a scale in a fixed metric and changes only by a bounded amount under a bi-Lipschitz change of that metric.

### 3.5 Distance entropy (Dai–Jiang 2008): a different notion with the same name

The literature already uses the phrase **distance entropy** for a different construction. For a separable metric space $X$, a continuous $f:X\to X$, a compact $K\subseteq X$ and $\epsilon>0$, the *step length* $n(f,K,\epsilon)$ is the largest $n$ such that $\operatorname{diam}(f^j(K))<\epsilon$ for all $0\le j<n$; set $\operatorname{diam}_\epsilon(K)=e^{-n(f,K,\epsilon)}$. Then a Hausdorff-type measure is built from covers by sets with prescribed step length, and its critical exponent is $\operatorname{ent}_H(f;Y)$, the distance entropy of $f$ on $Y$. Dai and Jiang prove that if $f$ is locally Lipschitz with constant $l$, then $\operatorname{ent}_H(f;Y)\le\max\{0,\mathrm{HD}(Y)\log l\}$, and if $f$ is locally expanding with skewness $\lambda$, then $\operatorname{ent}_H(f;Y)\ge\mathrm{HD}(Y)\log\lambda$; they also prove countable stability.

This is a **Hausdorff-dimension** construction using the time a set stays small, not the growth rate of a pair's distance. It is a valid "other choice" of how to turn distances into an entropy, but it is not your $h_{\mathrm{dist}}$. The name collision is worth flagging in the draft if you use "distance entropy" for $h_{\mathrm{dist}}$. Dai–Jiang's invariant is metric-bound in the same way — Hausdorff dimension and step lengths in a fixed metric — so it is not a topological invariant of $f$.

### 3.6 Length and volume growth (Newhouse, Yomdin)

The phrase "length growth rather than volume growth" in your abstract has a precise smooth analogue. For a $C^1$ map $f$ of a compact Riemannian manifold and $0\le k\le\dim M$, let
$$
v_k(f)=\sup_{D}\limsup_{n\to\infty}\frac1n\log\operatorname{vol}_k(f^n(D)),
$$
the supremum over $k$-dimensional disks $D$. Newhouse proved $v_k(f)\ge h_{\mathrm{top}}(f)$ for $C^r$ maps with $r>1$; Yomdin proved $v_k(f)\le h_{\mathrm{top}}(f)+k\,R(f)$ for $C^r$ maps, where $R(f)=\lim_n\frac1n\log\sup_x\|df^n_x\|$, and hence $v_\infty(f)=h_{\mathrm{top}}(f)$ for $C^\infty$ maps. For $k=1$ this is the exponential growth rate of lengths of curves; for $k=\dim M$ it is the volume growth.

The point for your draft is that **distance growth is the $k=0$ end of a dimension-indexed family**, with entropy at the $k=\dim M$ end. Your $h_{\mathrm{dist}}(x,y)$ is a pointwise $k=0$ growth rate; your $h_{\mathrm{diam}}$ and $h_{\max}$ are still $k=0$-type quantities (diameters and pair distances); $h_{\mathrm{top}}$ is the $k=\dim M$ volume growth. The Newhouse–Yomdin family gives a natural interpolation, and the "distance vs volume" dichotomy is the two ends of it. If you want an "other choice" that is intermediate between your $h_{\max}$ and $h_{\mathrm{top}}$, the $k$-volume growth rates are the natural candidates. The whole family is smooth: $v_k$ and $R$ are built from the derivative and the Riemannian volume, they are independent of the choice of metric, and they are defined only for smooth maps of manifolds (Newhouse's inequality needs $C^{1+\epsilon}$, Yomdin's $C^k$). The comparison with your $h_{\mathrm{dist}}$, $h_{\mathrm{diam}}$, $h_{\max}$ is an analogy between a smooth family and a topological one, not an identity.

A related "volume" notion is the **expansion entropy** of Hunt and Ott: the exponential growth rate of the average of the maximal $d$-dimensional volume growth ratio along orbits, minus the rate at which volume leaves the restraining region. On an invariant set it is the sum of the positive Lyapunov exponents; for $C^\infty$ maps it equals topological entropy (by Kozlovski). It generalizes Sacksteder–Shub's entropy and is the volume/Lyapunov-sum analogue of your distance-based expansion.

### 3.7 Lyapunov exponents and the growth spectrum

Your $h_{\mathrm{dist}}(x,y)$ is the exponential growth rate of the distance between the orbits of $x$ and $y$. For a $C^1$ map, this is a **pairwise Lyapunov exponent**: the growth rate of $d(f^n x,f^n y)$ is governed by the largest Lyapunov exponent of the pair (generically the larger of the two exponents, unless the pair is trapped in a stable direction). The minimum over pairs of this rate is the smallest Lyapunov exponent; the maximum is the largest. Thus your
$$
h_{\min}\le h_{\mathrm{diam}}\le h_{\max}\le h_{\mathrm{top}}
$$
is the topological analogue of
$$
\lambda_{\min}\le \text{(diameter growth)}\le \lambda_{\max}\le \sum_{\lambda_i>0}\lambda_i
$$
for smooth systems, with the last term being the volume growth (Pesin/Ruelle). The middle term depends on which notion of "diameter" you use; for a linear map it is the largest singular value, for a conformal map it is the common expansion rate.

The **Lyapunov spectrum** literature studies the level sets of the pointwise exponents: for $\alpha$, the set of points with exponent $\alpha$, its Hausdorff dimension, and the multifractal spectrum. Your "growth spectrum" is the set of values of $h_{\mathrm{dist}}(x,y)$; the natural refinement is the function $\alpha\mapsto$ (size of the set of pairs with $h_{\mathrm{dist}}=\alpha$). I did not find a paper that studies this for the cover-chain distance in the topological setting, but the smooth analogues are the multifractal analysis of Lyapunov exponents (Barreira–Pesin–Schmeling; Pesin–Weiss; Iommi–Jordan for interval maps). This is a possible "other choice" if you want the spectrum to be a function rather than a set. This whole paragraph is smooth: the exponents need a derivative and are not topological invariants. What survives topologically is the sign — positivity of the exponents, and the existence of a metric in which all of them are nonzero for an expansive homeomorphism (Fathi's metric; Section 8.4).

### 3.8 Visual metrics and conformal dimension

Your visual graph construction is a Floyd-style version of the **visual metric** on the boundary of a Gromov hyperbolic space. The Gromov product $\langle x,y\rangle_{\mathcal U}$ is the largest $n$ such that $x,y$ are $\mathcal U^{(n)}$-connected, and the visual metric is $\lambda^{-\langle x,y\rangle}$ for a discount factor $\lambda$. The standard references are M. Gromov, "Hyperbolic groups," in *Essays in Group Theory*, MSRI Publ. 8, Springer 1987, 75–263; M. Bonk and O. Schramm, "Embeddings of Gromov hyperbolic spaces," *Geom. Funct. Anal.* **10** (2000), 266–306; S. Buyalo and V. Schroeder, *Elements of Asymptotic Geometry*, EMS 2007. Floyd's original construction is W. J. Floyd, "Group completions and limit sets of Kleinian groups," *Invent. Math.* **57** (1980), 205–218. A note on kind: everything in this subsection is metric or conformal. Visual metrics depend on the discount factor, and $\mathrm{ARCdim}$ is invariant under quasisymmetric changes of metric but not under topological conjugacy.

The closest dynamical use is Haïssinsky–Pilgrim's **coarse expanding conformal** (cxc) theory. For a topologically cxc map of degree $d$, they construct visual metrics $d_\epsilon$ in which small balls are sent homeomorphically and homothetically onto their images with expansion constant $e^\epsilon$, and the Hausdorff dimension is $q=\log d/\epsilon$. The **Ahlfors regular conformal dimension** $\mathrm{ARCdim}$ is the infimum of Hausdorff dimensions over the conformal gauge; in the cxc setting it corresponds to the supremum of the admissible expansion exponents $\epsilon$, i.e. to $\log d/\mathrm{ARCdim}$. Your $\lambda_{\exp}$ is a supremum over **all** compatible metrics, not just the conformal gauge, so $\lambda_{\exp}\ge e^{\epsilon_{\max}}$ in the cxc case; whether equality holds is a separate question (and not one to settle here).

The thermodynamic formalism gives the general version of the equality: for a conformal expanding repeller, the Hausdorff dimension $s$ of the limit set is the root of Bowen's equation $P(-s\log|f'|)=0$, where $P$ is the topological pressure (Bowen 1979; Ruelle 1982; Parry–Pollicott 1990). This is the standard "distance vs volume" machine in the conformal case: the derivative (distance expansion) enters as a potential, and the pressure (volume growth) determines the dimension.

The energy side is your own paper with Pilgrim: K. M. Pilgrim and D. P. Thurston, "Ahlfors-regular conformal dimension and energies of graph maps," arXiv:2112.09041, *Groups Geom. Dyn.* (2024/2025), DOI 10.4171/ggd/914. The asymptotic $q$-conformal energy $\overline E^q[\pi,\phi]$ is continuous and non-increasing in $q$; $\mathrm{ARCdim}$ lies in the interval where $\overline E^q=1$. In the draft's comb-models section, the asserted identity $h_{\min}=-\log\overline E^\infty[\pi,\phi]$ is the $q=\infty$ end of this family, i.e. the maximal expansion rate, while $\mathrm{ARCdim}$ is the $q$ where the energy crosses 1. The two are different invariants, and it is worth keeping them explicitly separate: the maximal expanding constant is the reciprocal of the $\infty$-energy / stretch factor, not the reciprocal of the conformal dimension. (This is the user's own work; I am flagging the distinction, not claiming a new result. The draft's bibliography key for Haïssinsky–Pilgrim is `HP12:AlgExpanding`.)

---

## 4. Theme 2: the maximal expanding constant

### 4.1 Existence and equivalent definitions

The standard notion is: $f$ is **expanding** if there are $\lambda>1$ and $\epsilon>0$ such that $d(fx,fy)\ge\lambda d(x,y)$ whenever $d(x,y)<\epsilon$. Equivalent formulations in the literature:

- **Positively expansive + expanding metric.** Reddy proved that a map of a compact metric space is positively expansive iff it expands small distances iff it increases small distances. So every positively expansive map has an expanding metric.
- **Ruelle expanding.** A continuous open surjection that expands distances locally. On compact manifolds, positive expansiveness implies openness (invariance of domain), so the notions agree there.
- **Distance-expanding.** There are $\lambda>1$ and $C>0$ such that $\rho(f^n x,f^n y)\ge C\lambda^n\rho(x,y)$ for all $n$ and all $x,y$ sufficiently close (Przytycki–Urbański, *Conformal Fractals*, Ch. 4). This is the form used in thermodynamic formalism. Of the three, **positive expansiveness is topological** — it does not depend on the compatible metric — while **distance-expanding is metric-bound**: it is a property of one fixed metric. Ruelle's formulation sits between them: the openness is topological, the local expansion is metric. "Expanding in some compatible metric" is therefore a topological property even though "expanding in this metric" is not (Reddy, Section 2.1).

### 4.2 Fixed metric: the expansion constant and Hirsch's rate

For a map $f$ of a compact metric space and a compact positively invariant set $K$, the **expansion constant** is
$$
EC(f,K)=\sup\{\rho\ge0:\ \forall \rho^*<\rho\ \exists\delta>0,\ d(x,y)\le\delta,\ x\in K\implies d(fx,fy)\ge\rho^*d(x,y)\}.
$$
For $C^1$ $f$, $EC(f,K)=\min_{x\in K}m(Df_x)$, where $m(Df_x)=\inf_{|v|=1}|Df_xv|$ is the smallest singular value. Hirsch's **logarithmic expansion rate** is
$$
E(f,K)=\sup_{n\ge0}\frac1n\log EC(f^n,K).
$$
This is the maximal asymptotic rate at which $f$ expands distances in the fixed Riemannian metric. For a linear expanding map of $\mathbb T^n$, $EC(f)=\min_i|\lambda_i|$ and $E(f)=\log\min_i|\lambda_i|$. For a conformal expanding map, $EC(f)=\min|f'|$ and $E(f)=\log\lambda$ for the expansion rate $\lambda$. (The definition is quoted in Schreiber 1997.) A note on kind: $EC(f,K)$ is metric-bound — it is a property of the chosen metric — while its asymptotic rate $E(f,K)$ is not (Section 4.3). Neither is a topological invariant.

### 4.3 Schreiber's theorem: the smooth identification

Schreiber's theorem is the exact smooth analogue of your $\log\lambda_{\exp}=h_{\min}$. Let $f$ be a $C^1$ endomorphism of a Riemannian manifold, $K$ compact, positively invariant, and $f$ immersive at $K$ ($\det Df_x\ne0$ for $x\in K$). Then
$$
E(f,K)=\inf_\mu E(f,\mu),
$$
where $E(f,\mu)$ is the smallest Lyapunov exponent of $(f,\mu)$ and the infimum is over ergodic $f$-invariant measures supported on $K$. Equivalently,
$$
E(f,K)=\inf_{x\in K}\liminf_{n\to\infty}\frac1n\log m(Df^n_x).
$$
Consequences: $E(f,K)$ is independent of the Riemannian metric (Lyapunov exponents are), and it is invariant under bi-Lipschitz conjugacy. For a $C^{1+\alpha}$ map, $E(f,K)<-\lambda$ iff there are exponentially attracted points with rate $>\lambda$; for a one-dimensional map, $E(f,K)<-\lambda$ implies a periodic orbit with Lyapunov exponent $<-\lambda$. The paper also answers Hirsch's questions about aperiodic torus flows (expansion rate zero) and conjugacies (Hölder exponents). A note on kind: both sides of Schreiber's identity are smooth quantities — they need a differentiable structure and they are independent of the Riemannian metric. Your identity $\log\lambda_{\exp}=h_{\min}$ has the same shape but lives in the topological column on both sides; the analogy is one of form, not an instance of a common theorem.

This is the theorem to cite when you want to say "the maximal expanding constant is the infimum of the Lyapunov exponents" in the smooth case. In your topological setting, $h_{\min}$ replaces the infimum of the Lyapunov exponents, and $\log\lambda_{\exp}=h_{\min}$ is the analogue.

A measure-theoretic analogue is Morales (2025, arXiv:2504.16456): the expansion exponent of a measurable map is the minimum of the expansion exponents over Borel probability measures, and a map expands small distances (in Reddy's sense) iff every Borel probability measure has positive expansion exponent.

### 4.4 Over all metrics: your $\lambda_{\exp}$ and its nearest relatives

Your definition
$$
\lambda_{\exp}(f)=\sup\{\lambda>1:\ \exists\text{ compatible }d,\ \epsilon>0\text{ with }d(fx,fy)\ge\lambda d(x,y)\text{ for }d(x,y)<\epsilon\}
$$
is a supremum over **all** compatible metrics. Being a supremum over all compatible metrics, $\lambda_{\exp}$ is a topological conjugacy invariant by construction — the same device that makes $\lambda_{\mathrm{Lip}}$ topological (Section 8.5). I did not find this exact definition in the literature. The nearest relatives below are mostly theorems about particular metrics; only their suprema over a conjugacy-invariant class of metrics are topological. The nearest relatives are:

- **Dovbysh's optimal Lyapunov metrics** (2006). For expansive homeomorphisms, he constructs compatible Lyapunov metrics with prescribed contraction/expansion constants, and for A$^\#$-homeomorphisms proves that the lower bounds for the constants are attained simultaneously. This is the two-sided analogue of "choose the best metric." The one-sided/non-invertible case is not covered, but the techniques (Lyapunov metrics, Ruelle's Lipschitz axiom) are the same circle.
- **Fujita–Kato–Matsumoto's expanding ratios** (2010). For Ruelle expanding maps, they construct metrics with **exact** expansion ratio $\lambda$ for every $1<\lambda<s$, and prove the bound $\log\lambda\le h(f)/\dim X$. Their "strictly expands small distances" is stronger than your inequality expansion, so their set of exact ratios is a subset of the $\lambda$'s allowed by your definition; the suprema may or may not agree. Their Proposition 4.1 and Remark 4.4 are the cleanest published form of the entropy–dimension–expansion sandwich.
- **Roth–Roth's HausLip constant** (2020). They optimize $\mathrm{HD}_d(X)\log^+\operatorname{Lip}_d(f)$ over metrics and prove it equals $h(f)$ for positively expansive maps. Their construction gives a metric with $\operatorname{Lip}_D(f^n)\le16\,\mathrm{SkewLoc}_D(f^n)$, i.e. a near-uniform expansion metric, but it optimizes the product with the Hausdorff dimension, not $\log\lambda$ alone.
- **Fathi's hyperbolic metric** (1989). For hyperbolic sets, the metric that makes the stable/unstable rates uniform; the expanding case gives the entropy–dimension inequality.
- **Schreiber's $E$** (1997). In the smooth case, the supremum over Riemannian metrics of the asymptotic expansion rate is $E$, because $E$ is metric-independent. The one-step expansion constant $EC(f,K)$ does depend on the metric, but its asymptotic rate does not.

**One inequality between $\lambda_{\exp}$ and $E$ is elementary** (this is your observation, and I have not seen it stated). Let $d$ be a Riemannian metric and suppose $f^n$ is $(\rho,\epsilon)$-expanding in $d$, that is $EC_d(f^n,K)\ge\rho$. Put $\lambda=\rho^{1/n}$ and average the metric over the cycle:
$$
d'(x,y)=\sum_{i=0}^{n-1}\lambda^{-i}d(f^ix,f^iy).
$$
Then $d'$ is a compatible metric on a neighbourhood of $K$, and for $d(x,y)$ small,
$$
d'(fx,fy)-\lambda\,d'(x,y)=\lambda^{1-n}d(f^nx,f^ny)-\lambda\,d(x,y)\ \ge\ 0,
$$
because $d(f^nx,f^ny)\ge\rho\,d(x,y)=\lambda^n d(x,y)$. So $f$ is $(\lambda,\epsilon')$-expanding in $d'$, whence $\lambda_{\exp}\ge\rho^{1/n}$; taking the supremum over $n$,
$$
\lambda_{\exp}\ \ge\ e^{E(f,K)}.
$$
The class inclusion alone gives only $\lambda_{\exp}\ge\sup_{d\ \mathrm{Riem}}\lambda_{\mathrm{loc}}(d)$; the averaging step is what turns "one step of $f^n$" into "one step of $f$" and produces the exponent $1/n$. The reverse inequality does not follow this way, and I did not find it stated; by your theorem it is equivalent to $h_{\min}\le E(f,K)$.

The draft's theorem $\log\lambda_{\exp}=h_{\min}$ is, as far as I found, not in the literature. The closest statement is Schreiber's $E=\inf_\mu E(f,\mu)$; the draft's $h_{\min}$ is the topological replacement for the infimum of Lyapunov exponents. That is a clean way to position the theorem: it extends a known smooth identity to expansive topological systems, with the infimum over measures replaced by a minimum over pairs of points (or an in-radius growth rate).

### 4.5 The entropy–dimension–expansion sandwich

A note on kind: every quantity in this subsection is computed in a fixed metric — $D_d(X)$, $\mathrm{HD}_d(X)$ and $\log\lambda$ all depend on $d$. The inequalities matter because they bound the topological $h_{\mathrm{top}}$ by metric data, which is a bridge between columns rather than an identity inside one.

For a metric $d$ in which $f$ expands small distances by at least $\lambda$, the lower half of Fujita–Kato–Matsumoto's Proposition 4.1 (which uses only the lower bound) gives
$$
D_d(X)\log\lambda\le h_{\mathrm{top}}(f).
$$
Since the topological dimension is a lower bound for the box dimension, $D_d(X)\ge\dim X$, hence
$$
\log\lambda\le\frac{h_{\mathrm{top}}(f)}{\dim X},\qquad
\lambda_{\exp}(f)\le e^{h_{\mathrm{top}}(f)/\dim X}.
$$
This is the universal upper bound on the maximal expanding constant. It is the quantitative form of your $h_{\min}\le h_{\mathrm{top}}$: $\log\lambda_{\exp}=h_{\min}\le h_{\mathrm{top}}$, and the gap $h_{\mathrm{top}}/\dim X-h_{\min}$ measures the anisotropy of the expansion.

Conversely, if $f$ is Lipschitz with constant $L$ in the metric $d$, the Dai–Zhou–Geng/Misiurewicz inequality gives
$$
h_{\mathrm{top}}(f)\le \mathrm{HD}_d(X)\log L.
$$
For a metric that is both $\lambda$-expanding and $L$-Lipschitz, the two inequalities combine:
$$
D_d(X)\log\lambda\;\le\;h_{\mathrm{top}}(f)\;\le\;\mathrm{HD}_d(X)\log L.
$$
(Here $D_d(X)\ge\mathrm{HD}_d(X)$, so the two bounds are consistent.) If the metric is **exactly** $\lambda$-expanding (both a lower and an upper bound $\lambda$), then $h_{\mathrm{top}}(f)=D_d(X)\log\lambda$, which is the "fractal metric" formula of Fujita–Kato–Matsumoto. If it is only $\lambda$-expanding, the lower bound is the relevant one for $\lambda_{\exp}$; the upper bound involves the Lipschitz constant, which may be much larger.

The equality case $\lambda_{\exp}=e^{h_{\mathrm{top}}/\dim X}$ is what one expects for **conformal** expanding maps on a manifold: all directions expand at the same rate, the Hausdorff dimension of the adapted metric is the topological dimension, and Bowen's equation $h_{\mathrm{top}}=\dim\cdot\log\lambda$ holds. For non-conformal maps, the best metric can only expand the slowest direction at the slowest rate, so $\lambda_{\exp}<e^{h_{\mathrm{top}}/\dim X}$ in general. The 2-cat torus in your draft is the simplest example: $\lambda_{\exp}=\sqrt5-1\approx1.236$, while $e^{h_{\mathrm{top}}/2}=2$.

### 4.6 Examples where the maximal constant is known

In every case the invariant being computed is the topological $\lambda_{\exp}$; what changes from one example to the next is the structure used to compute it — smooth, algebraic, or combinatorial.

- **Circle doubling / graph maps.** For $z\mapsto z^d$ on $S^1$, $\lambda_{\exp}=d=e^{h_{\mathrm{top}}}$; for a graph endomorphism with irreducible Markov matrix, the draft's theorem gives $\log\lambda_{\exp}=h_{\min}=h_{\mathrm{top}}=\log\lambda_{\mathrm{Lip}}$, so $\lambda_{\exp}$ is the spectral radius of the Markov matrix. This is the dimension-one case where the upper bound is attained. For a general $C^1$ expanding map of $S^1$ of degree $d$, $h_{\mathrm{top}}=\log d$ (Misiurewicz–Przytycki), so $\lambda_{\exp}\le d$; equality holds for the constant-slope model, while a map such as $x\mapsto 2x+\epsilon\sin(2\pi x)$ has derivative varying between $2-2\pi\epsilon$ and $2+2\pi\epsilon$, so its maximal expanding constant is generally strictly less than $2$. The entropy–degree formula is Misiurewicz–Przytycki (1977).
- **Piecewise monotone interval maps.** Parry's constant-slope theorem: every transitive piecewise monotone interval map of positive entropy $h$ is conjugate to a map of constant slope $\pm e^h$ (Parry 1966). Hence the maximal expanding constant is $e^{h_{\mathrm{top}}}$ for such a map. Bobok–Roth show more generally that for every piecewise monotone interval map and every $C^\infty$ interval map the infimum of Lipschitz constants in the conjugacy class equals $e^{h_{\mathrm{top}}}$, and that the inequality can be strict for countably piecewise monotone maps; in the topologically mixing Markov case they identify the infimum with the exponential of the Salama entropy of a reverse Markov chain.
- **Linear torus maps.** For $A\in GL(n,\mathbb Z)$ with all eigenvalues $|\lambda_i|>1$, the best expanding constant is $\min_i|\lambda_i|$ (the slowest eigendirection), so $h_{\min}=\log\min_i|\lambda_i|$, while $h_{\mathrm{top}}=\sum_i\log|\lambda_i|$. The 2-cat matrix $\begin{pmatrix}0&2\\2&2\end{pmatrix}$ has eigenvalues $1\pm\sqrt5$, so $\lambda_{\exp}=\sqrt5-1$ and $h_{\mathrm{top}}=\log4$. This is my own elementary check, not a citation; it matches the draft's example.
- **Conformal expanding maps.** For a hyperbolic rational map on its Julia set, the visual metric gives a metric with $\lambda=e^{h_{\mathrm{top}}/\mathrm{HD}(J)}$ in the conformal gauge, so $\lambda_{\exp}\ge e^{h_{\mathrm{top}}/\mathrm{HD}(J)}$; the Ahlfors regular conformal dimension is the infimum of the dimensions of the gauge, and the maximal expansion exponent in the conformal gauge is $\log d/\mathrm{ARCdim}$ for a cxc map of degree $d$. The draft's companion paper on hyperbolic quadratic polynomials is the place where this is computed; I have not tried to reconcile its "period of the attracting cycle" statement with the conformal formula, since the two may use different normalizations.
- **Graph virtual endomorphisms.** The draft asserts $h_{\min}=-\log\overline E^\infty[\pi,\phi]$, i.e. $\lambda_{\exp}=1/\overline E^\infty[\pi,\phi]$, the reciprocal of the asymptotic $\infty$-conformal energy / stretch factor. This is the combinatorial analogue of Parry's constant-slope theorem and is the natural place where the maximal expanding constant is exactly computable. The stretch-factor background is Bestvina–Handel (1992) — the draft's `BH92:TrainTracksFreeGp` — and Thurston's 1986 norm paper; the draft's own keys for this material are `Thurston16:RubberBands`, `Thurston20:Characterize` and `PT26:ConfDim`.

### 4.7 Attainment and optimal metrics

Your $\lambda_{\exp}$ is a supremum, and the draft notes that it is appropriate because $\lambda$ can be decreased locally near a point. Whether the supremum is attained is a separate question. The literature on attainment:

- Dovbysh (2006) proves attainment of the lower bounds for the contraction/expansion constants for A$^\#$-homeomorphisms, with one metric satisfying all the conditions.
- Fujita–Kato–Matsumoto (2010) construct exact expanding metrics for every $\lambda$ below a bound, but they do not claim attainment at the supremum; their metrics approach the bound from below.
- Roth–Roth (2020) construct a metric with $\operatorname{Lip}$ and local skew comparable up to a factor 16, which is enough for the HausLip equality but not necessarily optimal for $\lambda$.
- Fathi (1989) constructs a hyperbolic metric for a hyperbolic set, again with uniform rates but not necessarily the optimal ones.

A natural question (not answered here) is whether $\lambda_{\exp}$ is attained for, say, graph maps or expanding maps of the circle. The examples above suggest it is attained in dimension one, but the general case is not covered by the papers I found.

---

## 5. Other choices in the draft's appendices

### 5.1 Neighborhoods of the diagonal and uniform spaces

Your Appendix `alternatives.tex` replaces open covers by open neighborhoods $N$ of the diagonal, with the distance defined by relation composition $N^k$. This is the **uniform space** formulation. The chain distance associated to an entourage is the standard "uniformity chain"; the volume versions are the packing/covering numbers of the uniformity. The entropy side is the **uniform entropy** of Bowen: for a uniform space $(X,\mathcal U)$ and a uniformly continuous self-map, one defines $r_n(V,K,\psi)$ as the minimal cardinality of a set $F$ such that $\bigcup_{x\in F}\bigcap_{k<n}\psi^{-k}(V(\psi^k(x)))$ covers $K$, and $h_{\mathrm{unif}}(\psi)=\sup_V\sup_K\limsup_n\frac1n\log r_n(V,K,\psi)$. For compact spaces this coincides with topological entropy; for noncompact spaces it is a uniform invariant. For the uniform entropy see Kimura (1998) and the survey of Dikranjan–Sanchis–Virili (2012); Bowen's entropy for uniformly continuous maps is Bowen (1971).

The point for your appendix: the neighborhood-of-the-diagonal language is not just a convenient reformulation; it is the standard language of uniform dynamics, and the "volume" definitions in that language are the uniform entropy ones. The distance definitions (relation composition $N^k$) are the uniform chain distance. The two-cover in-radius has a uniform analogue: the smallest $r(n)$ such that $(N^{(n)})^{r(n)}\subseteq N$, which is exactly the uniform version of your $r_{\mathcal U}(n)$.

### 5.2 Metrics, $\epsilon$-nets, and $\epsilon$-chains

Your second appendix's metric formulation is the standard Bowen/Dinaburg setup: $d_n^f(x,y)=\max_{k<n}d(f^k x,f^k y)$, $\epsilon$-separated/spanning sets, and $\epsilon$-chains. The "volume" notions are $\epsilon$-packing, $\epsilon$-covering, and $\epsilon$-nets; the "distance" notion is the $\epsilon$-chain length. The chain-based entropy formula of Misiurewicz (1986) is the relevant prior work if you want to recover the volume entropy from chains. The chain recurrence literature (Conley theory; Richeson–Wiseman) uses the same $\epsilon$-chains for recurrence and mixing times, and relates the chain mixing time to entropy.

A caution: the $\epsilon$-chain distance $d(x,y;\epsilon)$ is not comparable to the metric distance $d(x,y)$ from above unless the metric is a length metric. Your appendix notes this; the same caveat appears in the coarse geometry literature, where $\epsilon$-chains define a large-scale metric only after taking $\epsilon\to0$ and waiting for connectivity.

At the large-scale end there is a separate notion of **coarse entropy** (W. Geller, M. Misiurewicz, D. Sawicki, "Coarse entropy of metric spaces," *Geom. Dedicata* **218** (2024), article 99; arXiv:2109.15172). It uses $\delta$-pseudoorbits and $R$-separated sets with $\delta,R\to\infty$, is a coarse invariant of the metric space, and takes only the values $0$ or $\infty$. It is a distance-based entropy, but at large rather than small scales, so it is not a substitute for your $h_{\min}$ on a compact space; it is a natural "other choice" if you ever want a noncompact analogue.

### 5.3 Connected components of preimages

Your Remark in `intro.tex` suggests replacing the iterated cover $\mathcal U^{(n)}$ by the connected components of $f^{-(n+1)}(U)$ for $U\in\mathcal U$, at least for locally connected $X$. This is exactly the Haïssinsky–Pilgrim cxc setup: for a topologically cxc map, the preimages of a cover by connected open sets have connected components that shrink, and the visual metric is built from the resulting tree of components. The equivalence of the two constructions (iterated covers vs. connected components of preimages) is discussed in HP09, and the local connectivity hypothesis is what makes the components well-behaved. This is a good "other choice" to cite: it is not just a technical variant, it is the standard way to build visual metrics for cxc maps.

### 5.4 limsup/liminf vs limit

Your $h_{\mathrm{dist}}(x,y;\mathcal U)$ is defined with a limit and you note that it may not exist; one can replace it by limsup/liminf. Sun's $h_L^\pm$ and $h_E^\pm$ are the standard examples of this choice; Dai–Jiang's distance entropy also uses critical exponents rather than a limit. The draft's $h_{\min}$ uses Fekete's lemma (supermultiplicativity of $r_{\mathcal U}(n)$) to get a limit, which is the right move for the minimum growth. For $h_{\mathrm{diam}}$ and $h_{\max}$, the LEO hypothesis is what makes the limits exist; this is analogous to the standard "submultiplicative if LEO" arguments in complex dynamics.

### 5.5 Diameter growth, maximum growth, and LEO

Your $h_{\mathrm{diam}}$ and $h_{\max}$ are the "maximal" end of the distance spectrum. The smooth analogue is the maximal Lyapunov exponent / diameter growth of the manifold; the volume analogue is the Newhouse–Yomdin $k$-volume growth. The LEO (locally eventually onto) hypothesis is standard in one-dimensional dynamics and in the theory of expanding maps; it ensures that the diameter growth is coarsely submultiplicative. I did not find a paper that studies exactly your $h_{\mathrm{diam}}$ for topological dynamical systems, but the smooth and conformal cases are covered by the Lyapunov/volume-growth literature.

### 5.6 Polynomial and generalized growth

If $h_{\mathrm{top}}=0$, the exponential growth rate is zero and one can measure polynomial or other growth instead. The relevant notions are **polynomial entropy** (M. Marco, "Polynomial entropies and integrable Hamiltonian systems," *Regul. Chaotic Dyn.* **18** (2013), 623–655) and **generalized entropy** for arbitrary orders of growth (J. Correa and E. R. Pujals, "Orders of growth and generalized entropy," arXiv:2003.01257, 2020). The definitions are in the same spirit as yours: replace the exponential rate by a different growth scale, and measure the separation of orbits at that scale. This is a natural "other choice" if you want the theory to say something in the zero-entropy case.

### 5.7 Cleaning covers

Your `intro.tex` asks whether the operation — take unions of ascending chains, then delete dominated sets, then restrict to non-dominated elements — has a general name. The composite still has none that I can find, neither for covers nor in lattice theory. Asked separately, though, **both halves are standard, and each has a name**: the chain half is ordinary order theory, the deletion half is the theory of antichains. Neither name is a name *for covers*, and the cover literature's own reduction of a cover is a different and stronger operation. So the answer to "does this have some general name in lattice theory?" is: the vocabulary for both halves exists, the vocabulary for the composite does not.

**The chain half: adding upper bounds of ascending chains.** On a family of subsets, "adjoin the union of every ascending chain of members" is the operation of closing the family under unions of chains. In that generality it is a well-studied property, with three equivalent descriptions.

*Order theory.* A poset in which every chain has a least upper bound is **chain-complete**, and Markowsky (1976) is the paper that develops the notion; it proves the non-obvious direction, that a chain-complete poset is a **dcpo** — every directed subset has a least upper bound — so "closed under unions of chains" and "closed under unions of directed subfamilies" are one condition, not two (for families of sets this is Mayer-Kalkschmidt–Steiner 1964). The same paper gives the universal form of your operation: the **chain completion** of a poset is the free chain-complete poset on it, and Markowsky says explicitly that it is not the MacNeille completion. Two remarks that matter for your setting. First, nothing needs completing: $\mathcal O(X)$ is already a complete lattice, so the operation is the plain one of closing a subfamily under the chain suprema that already exist. Second, a chain-complete *lattice* is already a complete lattice (Markowsky's introduction), so for a lattice the chain hypothesis is not a weak one: it already yields all suprema.

*Universal algebra.* For families closed under intersections the dictionary is exact. Nation's Theorem 3.1: for a closure operator $\Gamma$, the following are equivalent — $\Gamma$ is **algebraic** (finitary); the union of any up-directed family of $\Gamma$-closed sets is $\Gamma$-closed; the union of any *chain* of $\Gamma$-closed sets is $\Gamma$-closed. So the chain step is exactly what makes a closure system algebraic, and the adjective for a family closed under directed unions in that literature is **inductive** (Erné's term).

*Convexity.* A family closed under arbitrary intersections and under directed unions is exactly a **convexity** in van de Vel's sense, and the theory of convex structures (van de Vel 1993) is where "closed under directed unions" is an axiom rather than a lemma.

*For covers specifically* I found no name and no use of this closure. Two standard facts sit next to it. First, in your setting the chain step is nearly vacuous at the top: $X$ is compact if and only if the top element of the frame $\mathcal O(X)$ is a compact element (Picado–Pultr 2012), i.e. if and only if every directed family of open sets whose union is $X$ already contains $X$; so if $X\notin\mathcal U$, no chain-union of members is $X$. Second, the classical Zorn argument for covers does take a union along an ascending chain, but the chain is a chain of *subcollections*, not of members: Arens and Dugundji (1950) call a subcollection redundant if the remaining members still cover, note that the union of an ascending chain of redundant subcollections is redundant — this is where point-finiteness enters — and conclude with Zorn's lemma that every open cover of a metacompact space has a minimal cover refining it.

**The deletion half: deleting dominated sets.** Keeping only the members contained in no other member produces what combinatorics calls a **Sperner family**, an **antichain**, or a **clutter**. The map $\max$ ("take the maximal members") needs no hypothesis to be *defined*, and it is always monotone, idempotent and contractive, i.e. an interior operator. What needs a hypothesis is **uniqueness**, and you were right to ask for one: the correspondence the counting literature states — down-set $\leftrightarrow$ antichain of its maximal members — is a bijection only under a chain condition, and finiteness is the special case usually meant.

- *What holds without hypotheses.* For *finitely many* antichains, "unite, then delete the dominated members" is the join in the order "every member of $\alpha$ lies in some member of $\beta$": $\alpha\vee\beta=\max(\alpha\cup\beta)$, written $\lceil\cdot\rceil$ in the antichain literature (De Wulf–Doyen–Henzinger–Raskin 2006). No finiteness of the ambient poset is needed here, because a chain in a union of $k$ antichains has at most $k$ elements, so nothing can be lost.
- *What needs the hypothesis.* Every member of the family must lie below a *maximal* member of the family. Finiteness, and more generally the ascending chain condition, gives this. In general the hypothesis is exactly the one you named — every ascending chain in the family has an upper bound in the family, which is Zorn's hypothesis. Without it the maximal set carries no information: on $[0,1]$, the families $\{(0,1-1/n)\}_{n\ge2}$ and $\{(0,1-1/n)\}_{n\ge2}\cup\{(0,1/2)\}$ both have empty maximal set, and the family $\{\{0,\dots,n\}\}_{n\ge0}$ of finite subsets of $\mathbb N$ has no maximal member at all. For infinitely many antichains the same failure kills the join: the antichains $\{\{0,\dots,n\}\}$ have no upper bound at all among the antichains of the finite subsets of $\mathbb N$, since no finite set contains every $\{0,\dots,n\}$, and therefore no least one.
- *Kleitman's conditions are exactly the finite case.* His objects are the subsets of a finite $n$-element set $S_n$, and his ideal is a down-set ($B\in I$ and $A\subseteq B$ implies $A\in I$). Then "an ideal can be uniquely corresponded to an antichain (its maximal members)" holds — in a finite poset every member of a down-set lies below a maximal member *of that down-set* — and the count is the Dedekind number. When $P$ is the Boolean lattice on $n$ atoms, that lattice of antichains is the **free distributive lattice on $n$ generators**, whose elements are the monotone Boolean functions (Birkhoff 1967, p. 63; Kleitman 1969). So Kleitman's sentence is your step *with a hypothesis supplied*, and the hypothesis he uses is finiteness, not a chain condition. That is why the sentence does not transfer to covers of a compact space as it stands: a cover of $[0,1]$ can contain an infinite ascending chain with no top, and then its antichain has lost a member.

For covers, the standard reduction is stronger, and this is the part worth getting right in the paper. The literature's notion is not "no member contains another" but **minimal**, or **irreducible**: no proper subfamily covers, equivalently no member is contained in the union of the *others* (Arens–Dugundji 1950; a space is called irreducible when every open cover has an irreducible open refinement). Deleting dominated members is strictly weaker. On a three-point space the cover $\{\{1,2\},\{2,3\},\{1,3\}\}$ has no member containing another, yet any two of its members cover. Irreducibility implies the Sperner condition and not conversely, so your deletion step is the *antichain* (Sperner) reduction of a cover, not the *irreducible* (minimal-subcover) reduction: on an irreducible cover the antichain reduction does nothing, and on an antichain cover it need not be irreducible. Nor is your deletion step "deleting redundant members" in the sense of Arens–Dugundji: a member can be redundant without being contained in any single other member, which is the same three-point example.

**What this means for the draft.** (1) If you want to name the composite rather than each half, "the antichain of maximal members of the chain-closure" is accurate but not a term of art; there is nothing to cite, and one sentence saying so is the right amount. (2) The two steps are not independent conveniences: the chain step is what makes the deletion step legitimate, by supplying precisely the chain condition above. From $\mathcal U$ alone the maximal members need not cover — on $[0,1]$, $\mathcal U=\{(0,1-1/n)\}_{n\ge2}\cup\{[0,1/2],\{1\}\}$ has maximal members $[0,1/2]$ and $\{1\}$, which miss $(1/2,1)$ — whereas with the chain-closure taken in the idempotent sense, Zorn's lemma puts every member below a maximal chain-union, so the maximal chain-unions do cover. (That argument is my reading of your definition, not a theorem quoted from a source.) (3) One caution about the definition as written: condition (1) quantifies over $\mathcal U$, not over the chain-closure, and then the output need not be an antichain. On $[0,1]$ take $\mathcal U=\{(0,1/2)\}\cup\{(1/j,1)\}_{j\ge2}\cup\{\{0\},\{1/2\},\{1\}\}$; it is a cover, and both $(0,1/2)$ and $(0,1)=\bigcup_j(1/j,1)$ satisfy your two conditions, nested. Quantifying condition (1) over the chain-closure makes the output exactly the set of maximal chain-unions, which is an antichain and, with the closure taken in the idempotent sense, still covers $X$. (4) A last simplification you can state without proof: a *finite* cover is already chain-complete, so for finite covers your cleaning is just the deletion of dominated sets.

---

## 6. What I did not find

To be explicit about the boundary of this search:

- I did **not** find a paper that defines your $\lambda_{\exp}$ as the supremum over **all compatible metrics** of the local expansion rate. The closest are Dovbysh's optimal Lyapunov metrics (two-sided, expansive homeomorphisms), Fujita–Kato–Matsumoto's exact expanding ratios (Ruelle expanding maps), and Schreiber's metric-independent asymptotic rate (smooth maps).
- I did **not** find a paper that proves your identity $\log\lambda_{\exp}=h_{\min}$ in the topological setting. The closest is Schreiber's $E=\inf_\mu E(f,\mu)$ for $C^1$ maps. This is the natural "prior work" boundary for the theorem.
- I did **not** find a paper that studies your pairwise growth spectrum $h_{\mathrm{dist}}(x,y)$ as a function of $(x,y)$ in the cover-chain setting. The smooth analogue is the Lyapunov spectrum / multifractal analysis; the cover analogue may be new.
- I did **not** find a paper that uses the in-radius $r(\mathcal U;\mathcal V)$ as a primary invariant. Sun's Lebesgue number is the closest, and the comparison $r\asymp\delta/\operatorname{mesh}$ is my own derivation, not a quoted theorem.
- I did **not** find a paper that defines your $h_{\mathrm{diam}}$ or $h_{\max}$ in the cover-chain sense, nor one that asserts or denies $h_{\max}=\log\lambda_{\mathrm{Lip}}$. The closest are Kifer's pointwise maximum/minimum characteristic exponents (metric-dependent, no supremum over metrics) and the adapted-metric theorems behind the smooth analogue; see §8.
- I did **not** find a paper that studies the infimum of the Lipschitz constant over *length* metrics, as opposed to over all compatible metrics (where it is trivially $\le1$), over a conjugacy class of Euclidean metrics (Bobok–Roth), or over metrics in a product with Hausdorff dimension (Roth–Roth).
- I did **not** find any paper that treats the composite operation of your `intro.tex` — close a cover under unions of ascending chains, then keep the maximal members — as an object with a name, or that uses it. The two halves are separately standard and separately named (Section 5.7); the composite is not.
- I did **not** read every paper in full. The statements above are from the abstracts, introductions, and the specific theorems I cite; I have not verified every proof. The Fujita–Kato–Matsumoto inequality is from their Proposition 4.1/Remark 4.4; the Sun inequalities are from the two Sun papers; the Schreiber theorem is from Schreiber's paper. The comparisons labelled "my own" are just that.

---

## 7. Follow-up: Julia sets, polynomial Julia sets, hyperbolic quadratic polynomials

This section answers the follow-up question specifically. It deliberately leaves out the conformal-dimension story (Haïssinsky–Pilgrim, ARCdim, the Pilgrim–Thurston energy paper), which is already known. I have separated work from your own group (Braga–Sieg, §7.1) from independent prior work (Hruska, §7.2; Gelfert–Przytycki–Rams, §7.3). A note on kind: everything in this section is smooth — metrics on $J$, derivatives, Lyapunov exponents. The one exception is the Braga–Sieg $\lambda(f)$: it is a supremum over metrics on $J$, so as an invariant of $(J,f)$ it is topological, even though every metric used to compute it is conformal.

### 7.1 Braga–Sieg: the exact sup-over-metrics definition on a Julia set (work from your group)

The closest match to your $\lambda_{\exp}$ on a Julia set is **Vitor Braga and David Sieg, “Expansivity on Julia sets”** (2023 Indiana University REU report, 20 pp., `https://math.indiana.edu/documents/reu/2023/Braga-Sieg.pdf`; presented at YMC 2024 as David R. Sieg, “Expansivity on Julia Sets”; I did not find a journal version). They define, for a rational map $f:\widehat{\mathbb C}\to\widehat{\mathbb C}$ with Julia set $J$:

- $f$ is $(\lambda,\epsilon)$-expansive on $J$ with respect to a metric $d$ if for all $z_1,z_2\in J$ with $d(z_1,z_2)<\epsilon$, one has $d(fz_1,fz_2)\ge \lambda d(z_1,z_2)$;
- $f$ is $\lambda$-expansive on $J$ if this holds for some $\epsilon>0$;
- the **expansivity** of $f$ on $J$ is
  $$
  \lambda(f)=\sup\{\lambda:\ \exists\text{ a metric }d\text{ on }J\text{ with }f\ (\lambda,\epsilon)\text{-expansive on }J\}.
  $$

This is exactly your $\lambda_{\exp}(f)$ restricted from $X$ to the Julia set $J$, with the same “sup over compatible metrics” and the same local inequality. It is the only place I found where this exact definition is written down for rational maps.

They relate $\lambda(f)$ to the asymptotic $\infty$-energy $E^\infty$ of the associated virtual endomorphism / spine. Their two key lemmas (stated in the report, with proofs deferred to a forthcoming research paper) are:

- **Lemma 14:** if $f^n$ has an $f^n$-invariant graph $G_0\subset G_n$ with induced map $g:G_0\to G_0$ and Markov partition matrix $A$, then
  $$
  \lambda(f)\ \le\ \bigl(\lambda_n(g)\bigr)^{1/n},
  $$
  the $n$-th root of the Perron–Frobenius eigenvalue.
- **Lemma 15:** if $E^\infty[\pi,\phi]<1$, then
  $$
  E^\infty[\pi,\phi]\ \ge\ \frac{1}{\lambda(f)}.
  $$
  With the invariant graph of Lemma 14, $E^\infty[\pi,\phi]\ge 1/(\lambda_n(g))^{1/n}$.

Thus invariant graphs give upper bounds for $\lambda(f)$, and $E^\infty$ gives lower bounds; when the two match, $\lambda(f)=1/E^\infty$. This is the same relation as the draft's asserted $h_{\min}=-\log E^\infty$ in the comb-models section, with $\log\lambda(f)=h_{\min}$.

**Quadratic examples.** They compute $E^\infty$ for several quadratic census maps of Brezin–Byrne–Levy–Pilgrim–Plummer:

- **Quadratic map #7:** $E^\infty=2^{-1/2}$, so the matching bounds give $\lambda(f)=\sqrt2$. (The report lists its critical portrait as type C(4,4); I have not decoded the notation, so I will not label the cycle period.)
- **Quadratic map #3.2** (period-3 critical cycle, the rabbit side): $E^\infty=2^{-1/3}$, so the matching bounds give $\lambda(f)=2^{1/3}$. This is exactly the “$\lambda_{\exp}$ records the period of the attracting cycle” phenomenon in a specific period-3 case.
- **Quadratic map #6.1** (a more complicated critical portrait): $E^\infty\approx0.739$, so the matching bounds give $\lambda(f)\approx1.353$ (if the displayed bounds are tight; the report's exact algebraic equation for this case is at the end of §4).
- **Cubic map #51:** the report only gets $0.553<E^\infty<0.563$; this is an example where the period alone does not determine the answer. (Lemma 15 turns the upper bound into $\lambda(f)>1.78$; an upper bound for $\lambda(f)$ would need a separate invariant-graph estimate.)

Caveats: this is an REU report, not a refereed paper; Lemmas 14–15 are announced with proofs deferred; and the examples are quadratic *census* maps, which include polynomial cases (the period-2 and period-3 critical-cycle cases above) but are stated for quadratic rational maps. Since Braga and Sieg were working under your supervision, this is your group's own work, not independent prior work. I record it here as the exact internal match to your definition; the independent comparison is Hruska in §7.2.

### 7.2 Hruska: explicit expanding metrics and box-expansion constants

The other direct line is **Suzanne Lynch Hruska, “Constructing an expanding metric for dynamical systems in one complex variable,”** *Nonlinearity* **18** (2005), no. 1, 81–100; arXiv:math/0406003; and her Cornell thesis, **“On the numerical construction of hyperbolic structures for complex dynamical systems”** (August 2002), available at `https://www.math.stonybrook.edu/theses/thesis02-2/part1.pdf`.

The paper's setup is the standard one: a polynomial $f:\mathbb C\to\mathbb C$ of degree $>1$ is hyperbolic if there is a Riemannian metric $\mu$ on a neighborhood $B$ of $J$ and an expansion constant $\lambda>1$ such that $|D_zf(v)|_\mu\ge\lambda|v|_\mu$ for every $z\in B$ and $v\in T_z\mathbb C$. Hruska introduces a **box metric**: a piecewise constant conformal metric $\Phi$ on a finite box cover $\{B_k\}$ of $J$, with “handicaps” $\phi_k>0$. The map is **box-expansive** by $L>1$ if, for every edge $B_k\to B_j$ of the box graph, every $z\in B_k$ and every $v\in T_z\mathbb C$,
$$
\phi_j\,|D_zf(v)|\ \ge\ L\,\phi_k\,|v|.
$$
Theorem 1.4 proves that box expansion by $L>1$ implies a smooth Riemannian metric on a neighborhood of $J$ with expansion at least $L$. So every computed $L$ is a **lower bound** for the maximal expanding constant of $f$ on $J$.

The key combinatorial result is Proposition 1.7: for a box graph $\Gamma$, let $\underline L$ be the minimum, over all simple cycles $B_0\to\cdots\to B_{n-1}\to B_0$, of the average cycle multiplier
$$
\bigl(\lambda_0\cdots\lambda_{n-1}\bigr)^{1/n},\qquad
\lambda_k=\min_{z\in B_k}|f'(z)|.
$$
Then a box metric with box expansion $L$ exists **iff** $L\le\underline L$. So $\underline L$ is the ideal box-expansion constant on that finite graph; as the grid refines, the minimum average cycle multiplier should converge to the infimum of Lyapunov exponents over periodic orbits, i.e. the smooth analogue of $e^{h_{\min}}$. I have not found a paper that states this convergence for Julia sets explicitly; it is the symbolic-image / Morse-spectrum convergence of Osipenko (Hruska cites Osipenko [17,18] as the related computational method).

**Quadratic examples.** Hruska computes explicit box-expansion constants (all lower bounds for $\lambda_{\exp}$ on $J$):

| Polynomial | Parameter | Box expansion constant $L$ | Source |
|---|---|---|---|
| basilica | $z^2-1$ | $1.14067$ | paper and thesis |
| aeroplane | $z^2-1.755$ | $1.05069$ | paper and thesis |
| Cantor Julia set | $z^2+0.35$ | $1.00778$ | paper |
| rabbit | $z^2-0.123+0.745i$ | $1.00508$ (subdivided $2^8$ grid) | thesis |
| double basilica | $z^2-1.31$ | $1.06758$ | thesis |
| five-eared rabbit | $z^2-0.504+0.563i$ | $1.0264$ | thesis |
| six-eared rabbit | $z^2+0.388+0.22i$ | not completed (memory limit) | thesis |

The published paper's examples are the basilica, aeroplane, and Cantor set; the thesis has the rabbit, double basilica, and five/six-eared rabbits. The computed constants are lower bounds, not exact $\lambda(f)$; they depend on the box model. The thesis's rabbit value $1.00508$ is much smaller than the $2^{1/3}\approx1.260$ suggested by the Braga–Sieg period-3 computation, which is consistent with it being a coarse lower bound rather than the supremum.

### 7.3 The Lyapunov-spectrum line

Outside the polynomial setting, the same “sup over visual metrics” question appears for expanding Thurston maps. Bonk–Meyer define visual metrics on the sphere and their expansion factors; Yin proves that the supremum of the expansion factors over visual metrics is the combinatorial expansion factor, that it is attained for Lattès maps, and that it need not be attained for Lattès-type maps (Q. Yin, “Lattès maps and combinatorial expansion,” *Ergodic Theory Dynam. Systems* **36** (2016), 1307–1342; Bonk–Meyer, *Expanding Thurston Maps*, AMS Math. Surveys Monogr. 225, 2017). This is the same “maximize expansion over a natural class of metrics” idea, but over visual metrics rather than all compatible metrics.

For a hyperbolic rational map, the range of Lyapunov exponents on $J$ is a natural “smooth” version of the distance-growth spectrum. **K. Gelfert, F. Przytycki, M. Rams, “On the Lyapunov spectrum for rational maps,”** *Math. Ann.* **348** (2010), no. 4, 965–1004; arXiv:0809.3363, prove that for a rational map with no critical points in $J$, the Lyapunov spectrum is described by the pressure function and has a range $[\alpha^-,\alpha^+]$. The lower endpoint $\alpha^-$ is the infimum of the Lyapunov exponents on $J$, so $e^{\alpha^-}$ is the maximal expansion rate in the smooth conformal category (the analogue of $\lambda_{\exp}$ for a hyperbolic polynomial). The paper applies to hyperbolic quadratic polynomials but does not compute the endpoint for specific parameters.

Measure-specific computations:

- **S. Sridharan and A. R. Tiwari, “The dependence of Lyapunov exponents of polynomials on its coefficients,”** arXiv:1906.03864 (2019): computes Lyapunov exponents of hyperbolic quadratic and cubic polynomials with respect to various Bernoulli measures, using thermodynamic formalism. This is not the infimum over all measures, but it is a computational line for expansion rates on Julia sets.
- **“Lyapunov exponents of polynomials with respect to certain weighted Lyubich's measures,”** arXiv:1908.07711: another measure-specific computation for hyperbolic quadratic polynomials.
- **The lower Lyapunov exponent at the critical value** (Levin–Przytycki–Shen, arXiv:1304.6178) is a different quantity: it is negative iff the unicritical polynomial has an attracting cycle. It is not the Julia-set expansion constant.

### 7.4 What I did not find

- I did **not** find a published journal paper that defines $\lambda(f)$ as the supremum over metrics on $J$ and computes it for hyperbolic quadratic polynomials. The only exact-definition match is the Braga–Sieg REU report (unpublished), which is your group's own work (supervised by you), not independent prior work.
- I did **not** find a paper that states the general formula $\lambda(f)=2^{1/p}$ for a hyperbolic quadratic polynomial with an attracting cycle of period $p$. The Braga–Sieg values $2^{-1/2}$ and $2^{-1/3}$ are consistent with it, and their “$2^{-1/p}\le E^\infty$” lower bound for a critical cycle of type $C(2,p)$ is the mechanism, but the general statement is not in the sources I found. This is likely the gap the draft's companion paper fills.
- Hruska's constants are **lower bounds**, not exact $\lambda(f)$; the exact values would require the ideal box constant in the limit, which the paper does not compute.
- I did not find any independent (non-REU, non-Thurston-circle) paper that uses the sup-over-metrics expansion constant on Julia sets. Hruska's box-expansion constant is the closest independent computational notion, and it is a lower bound rather than the supremum.

---

## 8. Follow-up: the other end — diameter growth, maximum growth, and the minimal Lipschitz constant

Your memory is right: the maximum-end material is in `clippings.tex` and never reached `min-growth.tex`. As you put it, it was formerly part of the paper and was cut when the planned scope changed, which is why it now reads as an unintegrated sketch — it carries a TODO list rather than an argument. It sits directly after your note "Start on maximum growth stuff, assuming LEO. This turns out not to be sufficient for submultiplicativity, and it's not the focus of this paper", so the definitions there are what the comparison below is against.

The short answer: there is a well-developed literature for the **pointwise** maximum exponent in a metric space (Kifer 1983 and its descendants), a small one for **infima of Lipschitz constants** (Bobok–Roth 2019, Roth–Roth 2020), and a large one for the **volume**-growth version (Newhouse, Yomdin). I found no paper that defines $h_{\mathrm{diam}}$ or $h_{\max}$ in the cover-chain sense, and none that asserts or denies $h_{\max}=\log\lambda_{\mathrm{Lip}}$.

### 8.1 What the draft has (quoted from `clippings.tex`)

- **Diameter growth.** "$h_{\mathrm{diam}}(f;\mathcal U)=\lim_{n\to\infty}\frac{\operatorname{diam}(X;\mathcal U_n)}{n}$ (if it exists)", with the note "(Something about refinement)". The displayed formula is missing the logarithm: if $\operatorname{diam}(X;\mathcal U_n)$ is the diameter in the chain distance $d(\cdot,\cdot;\mathcal U_n)$, that diameter grows like $\lambda^n$ in every expanding example, so the quotient without $\log$ diverges and its limit is not the growth rate; it reads as a typo for $\frac1n\log\operatorname{diam}(X;\mathcal U_n)$.
- **Maximum growth.** Recorded as a sketch rather than a definition: "Define $h_{\max}$ as supremum over covers $\mathcal U$ of diameter of $\mathcal U_n\cap U$, for $U\in\mathcal U_n$. Only finite if $X$ is locally connected."
- **Lipschitz stretch.** "$\lambda_{\mathrm{Lip}}$ is the infimum, over all length metrics on $X$, of the best Lipschitz constant for $f$", together with your remark that the length-metric restriction is necessary because snowflaking $d\mapsto d^\epsilon$ ($\epsilon<1$) decreases the Lipschitz constant.
- **The one inequality.** "$h_{\max}\le\log\lambda_{\mathrm{Lip}}$ when defined."
- **The chain.** $h_{\min}\le h_{\mathrm{diam}}\le h_{\max}\le h_{\mathrm{top}}$ for $f$ locally eventually onto, expansive, $X$ compact locally connected, all inequalities strict in examples; and the figure-eight example ($z\mapsto z^3$ on one loop, $z\mapsto z^2$ on the other) in which $\lambda_{\mathrm{Lip}}$ and $e^{h_{\max}}$ drop from $3$ to $2$ when passing to the $z\mapsto z^2$ subsystem.

### 8.2 Where the diameter analysis went

The comparison of $h_{\mathrm{diam}}$ with $h_{\max}$ — global versus localised growth, what the two give for a linear expanding endomorphism, and the 2Cat example — is analysis rather than prior work, so it is a separate note, `reports/diameter-growth-2cat.md`. Its conclusion, after several corrections: $h_{\mathrm{diam}}=\log2$ for the 2Cat with the ball cover (the map factors as $g\,h$ with $g=\operatorname{diag}(2,2)$ and $h\in GL_2(\mathbb Z)$; the chain-adjacency coincide exactly under $h^{n-1}$), while $h_{\max}=\log2\tau$ is the local eigen-expansion — so the four-value chain $h_{\min}<h_{\mathrm{diam}}<h_{\max}<h_{\mathrm{top}}$ of `clippings.tex` is intact and strict. What belongs in this memo is the literature, which is §8.3–8.7 below.

### 8.3 Kifer 1983: the closest published match

**Yuri Kifer, "Characteristic exponents of dynamical systems in metric spaces,"** *Ergodic Theory Dynam. Systems* **3** (1983), no. 1, 119–127; DOI 10.1017/S0143385700001838. For a one-parameter group $S^t$ of homeomorphisms of a separable metric space $X$ without isolated points (discrete or continuous time), he sets
$$
B_x(\delta,T)=\{y\in X\setminus\{x\}:\ d(S^tx,S^ty)<\delta\ \text{for all } t\in[0,T]\},
$$
$$
A_\delta(x,t)=\sup_{y\in B_x(\delta,t)}\frac{d(S^tx,S^ty)}{d(x,y)},\qquad
a_\delta(x,t)=\inf_{y\in B_x(\delta,t)}\frac{d(S^tx,S^ty)}{d(x,y)},
$$
and then
$$
\Lambda^\pm_\delta(x)=\limsup_{t\to\pm\infty}\frac{1}{|t|}\ln A_\delta(x,t),\qquad
\lambda^\pm_\delta(x)=\limsup_{t\to\pm\infty}\frac{1}{|t|}\ln a_\delta(x,t).
$$
Under an integrability hypothesis on a Borel invariant probability measure $\mu$ (satisfied whenever $S^t$ is uniformly Lipschitz for $|t|\le1$), Theorem 1 gives that $\mu$-almost everywhere the limsups are limits and that $\Lambda^+_\delta(x)=-\lambda^-_\delta(x)$, $\lambda^+_\delta(x)=-\Lambda^-_\delta(x)$; and if the hypothesis holds for all small $\delta$, the limits
$$
\Lambda^\pm(x)=\lim_{\delta\to0}\Lambda^\pm_\delta(x),\qquad \lambda^\pm(x)=\lim_{\delta\to0}\lambda^\pm_\delta(x)
$$
exist. Theorem 1 also gives the invariance of $\Lambda^\pm,\lambda^\pm$ under $S^t$ and their constancy $\mu$-a.e. when $\mu$ is ergodic. Theorem 2 does the same for an invariant set $Q$, with $B_Q(\delta,T)=\{y\in X\setminus Q: d(S^ty,Q)<\delta\ \forall t\in[0,T]\}$ and $d(y,Q)$ in the denominator — the exact analogue of localizing at a point, one level up. Theorem 3 is a stability criterion: $\Lambda^+(x)<0$ plus forward Lyapunov stability implies asymptotic stability, while $\Lambda^+(x)>0$ excludes forward stability. Theorem 4 and its Corollary identify the quantities in the smooth case:
$$
\Lambda^+(x)=\lim_{t\to\infty}\frac1t\ln\|DS^t_x\|,\qquad
\lambda^+(x)=\lim_{t\to\infty}\frac1t\ln\inf_{\xi\in T_x,\ \|\xi\|=1}\|DS^t_x\xi\|,
$$
$\mu$-a.e., so that "$\Lambda^+(x)$ and $\lambda^+(x)$ are the maximal and the minimal characteristic exponents at $x$". Remark 1 is the reason this is the right comparison for you: $\Lambda^+(x)>0$ and $\lambda^+(x)<0$ is "a very weak kind of hyperbolicity". Remark 2 matters for a one-sided setting: the positive-time half of Theorem 1 holds when $S^t$ is only a semigroup.

**How it matches the draft.** $A_\delta(x,t)$ is exactly "sup over $y$ of the distance ratio at time $t$, among the $y$ whose orbit shadows that of $x$ to accuracy $\delta$" — the metric counterpart of $\sup_y$ of the pairwise growth in $\mathcal U$-chains, with "shadows to accuracy $\delta$" replacing "lies in the same $\mathcal U_t$-atom". So $\sup_x\Lambda^+(x)$ is the metric-space form of $h_{\max}$, and $h_{\mathrm{diam}}$ is what one gets by moving the $\sup$ over $x$ outside (*and* dropping the shadowing restriction), which is a different and cruder object. The comparison between the chain and metric formulations is my reading, not a theorem in Kifer's paper — it is the same kind of comparison as $r(\mathcal U;\mathcal V)\asymp\delta(\mathcal U)/\operatorname{mesh}(\mathcal V)$ in §3, and it should be stated with the same hedge.

**What Kifer does not do.** He never takes a supremum over compatible metrics, so his $\Lambda^\pm$ are the analogues of your *fixed-cover* quantities, not of the topological invariants. And the quantity is genuinely metric-dependent: see §8.4. Kind: metric-bound — $\Lambda^\pm$ are unchanged by a bi-Lipschitz change of the metric and changed by a topological conjugacy; the topological content is the sign.

### 8.4 The pointwise-exponent line: metric dependence, and what survives it

**C. A. Morales, P. Thieullen, H. Villavicencio, "Lyapunov exponents on metric spaces,"** *Bull. Austral. Math. Soc.* **97** (2018), 153–162; DOI 10.1017/S0004972717000703. They replace Kifer's two limits by the pointwise Lipschitz constant: for $f:X\to X$ and $x$ not isolated,
$$
f'(x)=\limsup_{y\to x}\frac{d(fx,fy)}{d(x,y)},\qquad
\chi_f(x)=\limsup_{n\to\infty}\frac1n\log (f^n)'(x).
$$
Theorem 4.2: on a Riemannian manifold with the Riemannian distance, $\chi_f(x)=\Lambda_f(x)=\limsup_n\frac1n\log\|Df^n_x\|$ — the classical upper exponent. Theorem 3.2: at a Lyapunov-stable point $x$, $\chi_f(x)\le\Delta^+(f,x)$ (their rendering of Kifer's quantity); Remark 3.3 says they do not know whether the inequality can be strict. Theorem 5.1: a periodic point with $\chi_f(x)<0$ is asymptotically stable. Corollary 5.4: for an expansive homeomorphism of a compact metric space without isolated points there is a *compatible* metric in which $\inf_x\max\{\chi_f(x),\chi_{f^{-1}}(x)\}>0$ — Fathi's metric again.

Two corrections to the way this line is usually summarized, both of which matter for your definition:

- Their invariance is **bi-Lipschitz**, not topological. Theorem 4.4(3) says $\chi_{hfh^{-1}}=\chi_f\circ h^{-1}$ for their "diffeomorphisms", which are homeomorphisms $h$ with $\|h'\|_\infty<\infty$ and $\|(h^{-1})'\|_\infty<\infty$. Example 2.6 states the boundary explicitly: bi-Lipschitz equivalent metrics give the same $\chi_f$, and "this may be false if the metrics are only compatible". So $\chi_f$ is an invariant of the Lipschitz structure, not of the topological conjugacy class.
- The metric dependence is not a technicality. **M. J. Pacifico, J. L. Vieitez, "Lyapunov exponents for expansive homeomorphisms,"** *Proc. Edinburgh Math. Soc.* (2020), DOI 10.1017/S0013091519000579; arXiv:1704.05284, work in Kifer's framework and prove that on a compact Peano space an expansive homeomorphism has $\Lambda^+(x)\ge\gamma>0$ and $\lambda^+(x)\le-\gamma$ for some $\gamma>0$ and *every* $x$, in a suitable metric — while giving an example where the original metric does not see this. What expansiveness pins down is the **sign**, in some metric; not the value, and not in a fixed metric.

This is the precise sense in which a *topological* maximum-end invariant must either be built from covers (your $h_{\max}$) or be an infimum/supremum over metrics (your $\lambda_{\mathrm{Lip}}$). The pointwise exponents in the literature are neither, which is why they read as "almost" your objects.

### 8.5 $\lambda_{\mathrm{Lip}}$: why the length-metric restriction is doing real work, and the closest published relatives

**An elementary remark on the definition.** Without a restriction on the metric the definition collapses. For any continuous self-map of a compact metric space, any compatible $d$ and any $\epsilon>0$,
$$
d_\epsilon(x,y)=\sum_{m\ge0}(1+\epsilon)^{-m}d(f^mx,f^my)
$$
is a compatible metric and $\operatorname{Lip}_{d_\epsilon}(f)\le 1+\epsilon$. (Check: $d_\epsilon(fx,fy)=\sum_{m\ge0}(1+\epsilon)^{-m}d(f^{m+1}x,f^{m+1}y)=(1+\epsilon)\sum_{k\ge1}(1+\epsilon)^{-k}d(f^kx,f^ky)\le(1+\epsilon)d_\epsilon(x,y)$.) So the infimum over *all* compatible metrics is at most $1$ — for any map whatsoever — and it is $0$ for a constant map. The length-metric restriction is therefore not a technical convenience: it is the only thing preventing $\lambda_{\mathrm{Lip}}$ from being trivial, and your snowflake remark is a second face of the same fact. For a homeomorphism the same computation with weights $2^{-|n|}$ gives an invariant metric with Lipschitz constant $1$, again non-length in general. (This paragraph is my own computation, not a quotation.)

**Jozef Bobok, Samuel Roth, "The infimum of Lipschitz constants in the conjugacy class of an interval map,"** *Proc. Amer. Math. Soc.* **147** (2019), 255–269; arXiv:1612.06095. They study
$$
\Lambda(f)=\inf\{\lambda:\ \exists\ \text{a homeomorphism }\psi:[0,1]\to[0,1]\text{ with }\psi f\psi^{-1}\ \lambda\text{-Lipschitz}\},
$$
an infimum over the pullbacks of the Euclidean metric by homeomorphisms — not over length metrics, and not over all compatible metrics. Kind: $\Lambda(f)$ is a **topological** invariant, because the infimum runs over a conjugacy-invariant class of metrics; unlike a Lipschitz constant in a fixed metric, it sees only the conjugacy class. Your $\lambda_{\mathrm{Lip}}$ is topological for the same reason, with length metrics in place of pullbacks of the Euclidean metric. Results:

- $\Lambda(f)\ge e^{h(f)}$ for positive entropy (their §2; the general lower bound).
- $\Lambda(f)\le\nu(f)$ (Theorem 3.4), where $\nu(f)=\limsup_n(\operatorname{Var}f^n)^{1/n}$ is the growth rate of the **total variation** of the iterates — not the number of iterated preimages. (The number of preimages enters a different theorem; see the next point. I flag this because it is easy to conflate the two, and one summary of this paper does.)
- $\Lambda(f)=e^{h(f)}$ for piecewise monotone maps (Corollary 3.5) and for $C^\infty$ maps (Corollary 3.6).
- For countably piecewise monotone, topologically mixing, Markov maps,
$$
\log\Lambda(f)=h_{\mathrm{RevSal}}(\Sigma)=\limsup_{n\to\infty}\frac1n\log\#f^{-n}(\{x\})
$$
(Theorem 4.17), the Salama entropy of the reverse Markov chain, equal to the growth rate of the number of iterated preimages of a nearly arbitrary point. Example 4.19 shows this can be **strictly larger** than $h(f)$.
- If $f$ is locally eventually onto — their "leo" is your LEO for interval maps — then $\Lambda(f)=e^{h(f)}$ again (Theorem 4.15).

So the closest published relative of $\lambda_{\mathrm{Lip}}$ always satisfies $\ge e^{h}$, equals $e^{h}$ exactly when the combinatorics is tame (piecewise monotone, $C^\infty$, or leo), and exceeds $e^{h}$ only for maps with unbounded combinatorics. That is a suggestive parallel to your chain $h_{\min}\le h_{\mathrm{diam}}\le h_{\max}\le h_{\mathrm{top}}$, and it is only a parallel: the two infima are taken over different classes of metrics, so there is no implication either way.

**J. Milnor, *Dynamics: one-dimensional*, §13** (`https://www.math.stonybrook.edu/~jack/DYNOTES/dn13.pdf`). The one place I found where the "sup over metrics" move is stated as a programme, and it is stated at the *minimum* end: "Here the constant $c$ depends on the particular choice of Riemannian metric. However the expansion constant $k$ is independent of the metric. ... (To obtain a differentiable conjugacy invariant, we could take the supremum of all such expansion constants.)" It is a remark, not a theorem, but it is the smooth counterpart of your $\lambda_{\exp}$ and it points at the same design question at the other end.

**S. Roth, Z. Roth, "Inequalities for entropy, Hausdorff dimension, and Lipschitz constants,"** *Studia Math.* **250** (2020), 253–264; arXiv:1807.01196, already in §3, is the other infimum-over-metrics construction: $\mathrm{HausLip}(X,f)=\inf_d \mathrm{HD}_d(X)\log^+\operatorname{Lip}_d(f)$, a product rather than a bare Lipschitz constant, which is a topological conjugacy invariant and equals $h(f)$ for positively expansive maps and linear torus maps. Its method — construct a metric in which the Lipschitz constant is controlled uniformly over iterates — is the technology any proof of an inequality between $h_{\max}$ and $\log\lambda_{\mathrm{Lip}}$ would use.

### 8.6 The volume end: where the two ends meet in one inequality

**S. E. Newhouse, "Entropy and volume,"** *Ergodic Theory Dynam. Systems* **8** (1988), 283–299, and **Y. Yomdin, "Volume growth and entropy,"** *Israel J. Math.* **57** (1987), 285–300. Write $V_{\ell,k}(f)$ for the supremum of $\limsup_n\frac1n\log\omega(f^n\circ\sigma)$ over $C^k$ maps $\sigma$ of the unit $\ell$-cube into $M$, where $\omega$ counts $\ell$-volume with multiplicity, and $V_\ell(f)$ for the resulting $\ell$-volume growth rate; and set
$$
R(f)=\lim_{n\to\infty}\frac1n\log\max_{x\in M}\|Df^n(x)\|.
$$
Then (Newhouse) $h(f)\le V_{\ell}(f)$ for $C^{1+\epsilon}$ maps, and (Yomdin) $V_{\ell,k}(f)\le h(f)+\frac{2\ell}{k}R(f)$ for $C^k$ maps, so for $C^\infty$ maps the two together give $h(f)=V_\ell(f)$. Everything here is smooth: $V_{\ell,k}$ and $R$ are built from the derivative and the Riemannian volume, and are independent of the choice of Riemannian metric. (I quote the form in L.-S. Young's survey "Entropy in dynamical systems", in *Entropy*, Princeton Univ. Press, 2003, which states Newhouse's and Yomdin's theorems together.)

Two things are worth extracting. First, $R(f)$ **is** a diameter-growth quantity: $\frac1n\log$ of the maximum over $x$ of the operator norm at time $n$, i.e. a maximum over $x$ taken *at each $n$* and then a limit — a $\lim_n\sup_x$, not a $\sup_x\lim_n$. That is the same order of quantifiers as your $h_{\max}$ (a supremum over the atoms of $\mathcal U_n$ at each $n$, then a limit), and it is not the same as $\sup_x\Lambda^+(x)$; the two differ exactly when the limsup is not uniform in $x$. Second, Yomdin's inequality is the one place in the literature where entropy and a maximum-end quantity appear in the same inequality: $V_{\ell,k}\le h+\frac{2\ell}{k}R$. It says the volume growth is bounded by entropy plus the maximum distance growth, and for $C^\infty$ maps the volume end collapses onto entropy. That is the closest published analogue I found of your chain, though it runs in the other direction (volume above entropy, bounded by entropy plus the max end).

### 8.7 What is missing from the literature

- **No cover-chain definitions.** I did not find a paper defining $h_{\mathrm{diam}}$ or $h_{\max}$ as you do. The phrase "diameter growth" is unusable as a search term: it returns Abresch–Gromoll and Riemannian geometry, Schreier-graph diameter growth in group theory, and manifold topology. "Maximal Lyapunov exponent", "maximal characteristic exponent" and "diameter of iterates" are the working terms, and they lead to Kifer.
- **No comparison theorem.** Nothing I found states $h_{\max}=\log\lambda_{\mathrm{Lip}}$, in either direction beyond your Proposition. Kifer and Morales–Thieullen–Villavicencio give pointwise, metric-dependent exponents; Bobok–Roth give the infimum in a conjugacy class of interval maps; Roth–Roth give an infimum-over-metrics product. None of them optimizes over length metrics, and none compares the result to a cover-based maximum growth.
- **No comparison theorem, and no counterexample either.** Nothing found states $h_{\max}=\log\lambda_{\mathrm{Lip}}$ in either direction. What the missing inequality would be, and where a counterexample would have to live, are set out in `reports/diameter-growth-2cat.md`, §7.

---

## 9. Bibliography

Every work cited in this memo, collected once and grouped by thread. Within each theme the order is alphabetical by first author. A work that belongs to two threads is listed under the primary one and named in the other theme's lead sentence. Passing mentions in the text that carry no year — the multifractal papers of Barreira–Pesin–Schmeling, Pesin–Weiss and Iommi–Jordan in §3.7, Kozlovski's theorem on expansion entropy, Osipenko's computational method, Conley's chain recurrence — are not indexed here.

**A. Entropy from distances: the classical layer**

- **R. L. Adler, A. G. Konheim, M. H. McAndrew**, "Topological entropy," *Trans. Amer. Math. Soc.* **114** (1965), 309–319.
- **R. Bowen**, "Entropy for group endomorphisms and homogeneous spaces," *Trans. Amer. Math. Soc.* **153** (1971), 401–414.
- **R. Bowen**, "Topological entropy for noncompact sets," *Trans. Amer. Math. Soc.* **184** (1973), 125–136.
- **E. I. Dinaburg**, "On the relations among various entropy characteristics of dynamical systems," *Math. USSR Izv.* **5** (1971), 337–378.
- **W. H. Gottschalk, G. A. Hedlund**, *Topological Dynamics*, AMS Colloquium Publications **36**, 1955.
- **P. Walters**, *An Introduction to Ergodic Theory*, Graduate Texts in Mathematics **79**, Springer, 1982.

**B. Lebesgue numbers, expansive constants, in-radii**

- **F. Przytycki, M. Urbański**, *Conformal Fractals: Ergodic Theory Methods*, London Math. Soc. Lecture Note Series **371**, Cambridge Univ. Press, 2010 (Ch. 4, including §4.6, "Expansive maps are expanding in some metric").
- **W. L. Reddy**, "Expanding maps on compact metric spaces," *Topology Appl.* **13** (1982), 327–334.
- **M. Shub**, "Endomorphisms of compact differentiable manifolds," *Amer. J. Math.* **91** (1969), 175–199.
- **P. Sun**, "Exponential decay of Lebesgue numbers," *Discrete Contin. Dyn. Syst.* **32** (2012), 3773–3785; arXiv:1002.0549.
- **P. Sun**, "Exponential decay of expansive constants," arXiv:1101.0851 (2011).

**C. Chains, uniform spaces, and other growth scales**

- **J. Correa, E. R. Pujals**, "Orders of growth and generalized entropy," arXiv:2003.01257 (2020).
- **X. Dai, Y. Jiang**, "Distance entropy of dynamical systems on noncompact-phase spaces," *Discrete Contin. Dyn. Syst.* **20** (2008), 313–333; DOI 10.3934/dcds.2008.20.313.
- **D. Dikranjan, M. Sanchis, S. Virili**, "New and old facts about entropy in uniform spaces and topological groups," *Topology Appl.* **159** (2012), no. 7; DOI 10.1016/j.topol.2011.05.046.
- **W. Geller, M. Misiurewicz, D. Sawicki**, "Coarse entropy of metric spaces," *Geom. Dedicata* **218** (2024), art. 99; arXiv:2109.15172.
- **T. Kimura**, "Completion theorem for uniform entropy," *Comment. Math. Univ. Carolin.* **39** (1998).
- **M. Marco**, "Polynomial entropies and integrable Hamiltonian systems," *Regul. Chaotic Dyn.* **18** (2013), 623–655.
- **M. Misiurewicz**, "Remark on the definition of topological entropy," in *Dynamical Systems and Partial Differential Equations* (Caracas, 1984), Univ. Simón Bolívar, Caracas, 1986, 65–67.
- **D. Richeson, J. Wiseman**, "Chain recurrence rates and topological entropy," *Topology Appl.* **156** (2008), 251–261; arXiv:0801.1635.

**D. The maximal expanding constant**

- **M. Bestvina, M. Handel**, "Train tracks and automorphisms of free groups," *Ann. of Math.* **135** (1992), 1–51.
- **X. Dai, Z. Zhou, X. Geng**, "Some relations between Hausdorff-dimensions and entropies," *Sci. China Ser. A* **41** (1998), 1068–1075.
- **S. A. Dovbysh**, "Optimal Lyapunov metrics of expansive homeomorphisms," *Izv. Math.* **70** (2006), 883–929; DOI 10.1070/IM2006v070n05ABEH002332.
- **A. Fathi**, "Expansiveness, hyperbolicity and Hausdorff dimension," *Comm. Math. Phys.* **126** (1989), 249–262.
- **C. Fujita, H. Kato, M. Matsumoto**, "Fractal metrics of Ruelle expanding maps and expanding ratios," *Topology Appl.* **157** (2010), 615–628; DOI 10.1016/j.topol.2009.11.002.
- **M. W. Hirsch**, "Asymptotic phase, shadowing, and reaction-diffusion systems," in *Differential Equations, Dynamical Systems, and Control Science*, Lecture Notes in Pure and Appl. Math. **152**, Dekker, 1994, 87–99.
- **C. A. Morales**, "A measure-theoretic expansion exponent," arXiv:2504.16456 (2025).
- **M. Misiurewicz**, "On Bowen's definition of topological entropy," *Discrete Contin. Dyn. Syst.* **10** (2004), 827–833.
- **M. Misiurewicz, F. Przytycki**, "Topological entropy and degree of smooth mappings," *Bull. Acad. Polon. Sci. Sér. Sci. Math. Astronom. Phys.* **25** (1977), 573–574.
- **M. Misiurewicz, W. Szlenk**, "Entropy of piecewise monotone mappings," *Studia Math.* **67** (1980), 45–63.
- **W. Parry**, "Symbolic dynamics and transformations of the unit interval," *Trans. Amer. Math. Soc.* **122** (1966), 368–378.
- **S. Roth, Z. Roth**, "Inequalities for entropy, Hausdorff dimension, and Lipschitz constants," *Studia Math.* **250** (2020), 253–264; arXiv:1807.01196.
- **S. J. Schreiber**, "Expansion rates and Lyapunov exponents," *Discrete Contin. Dyn. Syst.* **3** (1997), 433–438; DOI 10.3934/dcds.1997.3.433.
- **W. P. Thurston**, "A norm for the homology of 3-manifolds," *Mem. Amer. Math. Soc.* **59** (1986), no. 339.

The definitions this theme rests on — Reddy (1982), Przytycki–Urbański (2010), Shub (1969) — are listed under B. The closest relative of $\lambda_{\mathrm{Lip}}$, Bobok–Roth (2019), is under E, and Fathi (1989) is also used in G.

**E. The maximum end: diameter growth, maximum growth, and the minimal Lipschitz constant**

- **J. Bobok, S. Roth**, "The infimum of Lipschitz constants in the conjugacy class of an interval map," *Proc. Amer. Math. Soc.* **147** (2019), 255–269; arXiv:1612.06095.
- **Y. Kifer**, "Characteristic exponents of dynamical systems in metric spaces," *Ergodic Theory Dynam. Systems* **3** (1983), no. 1, 119–127; DOI 10.1017/S0143385700001838.
- **J. Milnor**, *Dynamics: one-dimensional*, Stony Brook lecture notes, §13 ("Expanding maps"), `https://www.math.stonybrook.edu/~jack/DYNOTES/dn13.pdf`.
- **C. A. Morales, P. Thieullen, H. Villavicencio**, "Lyapunov exponents on metric spaces," *Bull. Austral. Math. Soc.* **97** (2018), 153–162; DOI 10.1017/S0004972717000703.
- **M. J. Pacifico, J. L. Vieitez**, "Lyapunov exponents for expansive homeomorphisms," *Proc. Edinburgh Math. Soc.* (2020); DOI 10.1017/S0013091519000579; arXiv:1704.05284.

The other infimum-over-metrics construction, Roth–Roth (2020), is listed under D.

**F. Volume growth and Lyapunov exponents**

- **B. R. Hunt, E. Ott**, "Defining chaos," *Chaos* **25** (2015), 097618; arXiv:1501.07896.
- **S. E. Newhouse**, "Entropy and volume," *Ergodic Theory Dynam. Systems* **8** (1988), 283–299.
- **V. I. Oseledets**, "A multiplicative ergodic theorem: Lyapunov characteristic numbers for dynamical systems," *Trans. Moscow Math. Soc.* **19** (1968), 197–221.
- **Ya. B. Pesin**, "Lyapunov characteristic exponents and smooth ergodic theory," *Russian Math. Surveys* **32** (1977), 55–114.
- **R. Sacksteder, M. Shub**, "Entropy of a differentiable map," *Adv. Math.* **28** (1978), 181–185.
- **Y. Yomdin**, "Volume growth and entropy," *Israel J. Math.* **57** (1987), 285–300.
- **L.-S. Young**, "Entropy in dynamical systems," in *Entropy*, Princeton Univ. Press, 2003.

**G. Visual metrics, conformal dimension, and energies**

- **M. Bonk, O. Schramm**, "Embeddings of Gromov hyperbolic spaces," *Geom. Funct. Anal.* **10** (2000), 266–306.
- **R. Bowen**, "Hausdorff dimension of quasicircles," *Publ. Math. Inst. Hautes Études Sci.* **50** (1979), 11–25.
- **S. Buyalo, V. Schroeder**, *Elements of Asymptotic Geometry*, EMS Monographs in Mathematics, European Mathematical Society, 2007.
- **W. J. Floyd**, "Group completions and limit sets of Kleinian groups," *Invent. Math.* **57** (1980), 205–218.
- **M. Gromov**, "Hyperbolic groups," in *Essays in Group Theory*, MSRI Publications **8**, Springer, 1987, 75–263.
- **P. Haïssinsky, K. M. Pilgrim**, "Coarse expanding conformal dynamics," *Astérisque* **325** (2009).
- **W. Parry, M. Pollicott**, "Zeta functions and the periodic orbit structure of hyperbolic dynamics," *Astérisque* **187–188** (1990).
- **K. M. Pilgrim, D. P. Thurston**, "Ahlfors-regular conformal dimension and energies of graph maps," arXiv:2112.09041; *Groups Geom. Dyn.*, DOI 10.4171/ggd/914.
- **D. Ruelle**, "Repellers for real analytic maps," *Ergodic Theory Dynam. Systems* **2** (1982), 99–107.

**H. Julia sets and expanding metrics in one complex variable**

- **M. Bonk, D. Meyer**, *Expanding Thurston Maps*, AMS Mathematical Surveys and Monographs **225**, 2017.
- **V. Braga, D. Sieg**, "Expansivity on Julia sets," 2023 Indiana University REU report, `https://math.indiana.edu/documents/reu/2023/Braga-Sieg.pdf`.
- **K. Gelfert, F. Przytycki, M. Rams**, "On the Lyapunov spectrum for rational maps," *Math. Ann.* **348** (2010), 965–1004; arXiv:0809.3363.
- **S. L. Hruska**, "Constructing an expanding metric for dynamical systems in one complex variable," *Nonlinearity* **18** (2005), 81–100; arXiv:math/0406003; and the thesis *On the numerical construction of hyperbolic structures for complex dynamical systems*, Cornell University, 2002, `https://www.math.stonybrook.edu/theses/thesis02-2/part1.pdf`.
- **G. Levin, F. Przytycki, W. Shen**, "The Lyapunov exponent of holomorphic maps," arXiv:1304.6178.
- **S. Sridharan, A. R. Tiwari**, "The dependence of Lyapunov exponents of polynomials on its coefficients," arXiv:1906.03864 (2019).
- **S. Sridharan, A. R. Tiwari**, "Lyapunov exponents of polynomials with respect to certain weighted Lyubich's measures," arXiv:1908.07711 (2019).
- **Q. Yin**, "Lattès maps and combinatorial expansion," *Ergodic Theory Dynam. Systems* **36** (2016), 1307–1342.

**I. Covers, refinement, and uniform spaces**
- **R. Engelking**, *General Topology*, revised and completed edition, Sigma Series in Pure Mathematics **6**, Heldermann Verlag, Berlin, 1989.
- **J. R. Isbell**, *Uniform Spaces*, Mathematical Surveys **12**, American Mathematical Society, Providence, RI, 1964.
- **J. L. Kelley**, *General Topology*, D. Van Nostrand, Princeton, NJ, 1955; reprinted as Graduate Texts in Mathematics **27**, Springer, 1975.
- **S. Willard**, *General Topology*, Addison-Wesley, Reading, MA, 1970.

**J. Metric entropy and the partition order**
- **A. N. Kolmogorov**, "A new metric invariant of transitive dynamical systems and automorphisms of Lebesgue spaces," *Dokl. Akad. Nauk SSSR* **119** (1958), 861–864.
- **V. A. Rokhlin**, "Lectures on the entropy theory of measure-preserving transformations," *Russian Math. Surveys* **22** (1967), no. 5, 1–52; Russian original in *Uspekhi Mat. Nauk* **22** (1967), no. 5, 3–56.
- **Ya. G. Sinai**, "On the concept of entropy of a dynamical system," *Dokl. Akad. Nauk SSSR* **124** (1959), 768–771.

**K. Covers as a lattice: chains, antichains, and the two halves of "cleaning"**

- **R. Arens, J. Dugundji**, "Remark on the concept of compactness," *Portugaliae Math.* **9** (1950), 141–143.
- **G. Birkhoff**, *Lattice Theory*, 3rd ed., American Mathematical Society Colloquium Publications **25**, Providence, RI, 1967 (free distributive lattices and antichains, p. 63).
- **M. De Wulf, L. Doyen, T. A. Henzinger, J.-F. Raskin**, "Antichains: a new algorithm for checking universality of finite automata," in *Computer Aided Verification (CAV 2006)*, Lecture Notes in Comput. Sci. **4144**, Springer, 2006, 17–30.
- **D. J. Kleitman**, "On Dedekind's problem: the number of monotone Boolean functions," *Proc. Amer. Math. Soc.* **21** (1969), no. 3, 677–682; DOI 10.1090/S0002-9939-1969-0241334-6; MR 0241334.
- **G. Markowsky**, "Chain-complete posets and directed sets with applications," *Algebra Universalis* **6** (1976), 53–68.
- **J. Mayer-Kalkschmidt, E. Steiner**, "Some theorems in set theory and applications in the ideal theory of partially ordered sets," *Duke Math. J.* **31** (1964), 287–289.
- **J. B. Nation**, *Notes on Lattice Theory*, Ch. 3 (Theorem 3.1), `https://math.hawaii.edu/~jb/math618/os3uh.pdf`.
- **J. Picado, A. Pultr**, *Frames and Locales: Topology without Points*, Frontiers in Mathematics, Birkhäuser/Springer, Basel, 2012.
- **M. L. J. van de Vel**, *Theory of Convex Structures*, North-Holland Mathematical Library **50**, North-Holland, Amsterdam, 1993.
