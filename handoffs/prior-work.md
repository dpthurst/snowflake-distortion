*Copied into this repository verbatim on 2026-11-01, from `handoffs/prior-work.md` in the Sofic presentations project. It describes that session and its paths, not this repository, and is superseded: the memo now lives here in `reports/`, with `refs/` and the build script alongside it.*

# Handoff — prior-work thread (the *Expansion* literature memo)

*Prepared by Wren, an Open Math OS Agent. For a session continuing the **literature** side of this work, not the diameter side (that is `handoffs/diameter.md`).*

## 1. What this thread is

Dylan Thurston (DPT) has an early draft, `dpthurst/expansion` (GitHub), about quantifying expansion in topological dynamical systems. Two threads in the draft:

1. **Distance-growth analogues of topological entropy** — the draft's $h_{\min}$, $h_{\mathrm{dist}}(x,y)$, in-radius $r(\mathcal U;\mathcal V)$.
2. **Quantifying the maximal expanding constant** — the draft's $\lambda_{\exp}$.

DPT's framing, verbatim: *"I don't want you to think about any of the open questions there, but I am interested in any related prior work"*, and *"Be aware that the terminology will probably differ."* This thread is **literature only** — no proving, no resolving the draft's conjectures.

## 2. What exists, and where

| File | What it is |
|---|---|
| `reports/expansion-prior-work.md` | **The deliverable.** Literature memo, 28 pp as PDF. Sections 1–9. |
| `reports/expansion-prior-work.tex` | Generated from the `.md`; do not hand-edit. |
| `reports/expansion-prior-work.pdf` | Compiled output, **untracked** (`.gitignore` has `reports/*.pdf`). |
| `scripts/build-report-tex.py` | Builds either report: `python3 scripts/build-report-tex.py [name]` → `.tex` + `.pdf` in one run. |

Memo structure: §1 six papers to read first; §2 objects in the literature's language, **§2.1 the three kinds of invariant** (topological / metric-bound / smooth) plus three terminology warnings; §3 Theme 1 (Bowen–Dinaburg, $\epsilon$-chains, Sun's Lebesgue numbers, Sun's expansive constants, Dai–Jiang, Newhouse–Yomdin, Lyapunov spectra, visual metrics); §4 Theme 2 (existence, Hirsch's $EC$ and $E$, Schreiber, $\lambda_{\exp}$ over all metrics, the entropy–dimension sandwich, examples, attainment); §5 other choices in the draft's appendices; §6 what was not found; §7 Julia-set follow-up; §8 the maximum end (Kifer, the pointwise line, $\lambda_{\mathrm{Lip}}$, the volume end, gaps); §9 bibliography (78 entries, themes A–K, every one cited in the body).

## 3. What is established (check before restating)

- **The closest published invariant to $h_{\min}$ is Sun's $h_L$** (Lebesgue-number decay), DCDS 32 (2012); and $r(\mathcal U;\mathcal V)\asymp\delta(\mathcal U)/\operatorname{mesh}(\mathcal V)$ is *my* derivation, not a quoted theorem. Sun's $h_L$ is **metric-bound** (bi-Lipschitz invariant), *not* topological — his Examples 6.4–6.5 are topologically conjugate homeomorphisms with different $h_L$.
- **Schreiber's theorem** $E(f,K)=\inf_\mu E(f,\mu)$ (DCDS 3, 1997) is the smooth analogue of $\log\lambda_{\exp}=h_{\min}$.
- **$\lambda_{\exp}\ge e^{E(f,K)}$** — DPT's observation, carried out in §4.4: if $f^n$ is $(\rho,\epsilon)$-expanding in a Riemannian $d$, averaging $d'(x,y)=\sum_{i<n}\rho^{-i/n}d(f^ix,f^iy)$ gives $d'(fx,fy)-\rho^{1/n}d'(x,y)=\rho^{(1-n)/n}d(f^nx,f^ny)-\rho^{1/n}d(x,y)\ge0$. The reverse inequality is open and, by the draft's theorem, equivalent to $h_{\min}\le E$.
- **"Cleaning" a cover (memo §5.7)**: the *composite* operation has no name in the literature, but each half does, and the memo now says so. The chain half — adjoining unions of ascending chains — is **chain-completeness**: Markowsky 1976 (chain-complete = dcpo; the *chain completion*, which is not the MacNeille completion), Mayer-Kalkschmidt–Steiner 1964 (chain unions = directed unions for families of sets), Nation Thm 3.1 (algebraic/finitary closure systems are exactly the chain-closed ones), van de Vel 1993 (intersections + directed unions = a convexity). The deletion half is the **antichain (Sperner, clutter) reduction**: the maximal-elements map is an interior operator and "unite then delete dominated" is the join in the lattice of antichains; on a Boolean lattice that lattice is the free distributive lattice (Birkhoff 1967; Kleitman 1969). For covers, the standard reduction is the *stronger* one: **irreducible / minimal** = no proper subcover, equivalently no member inside the union of the others (Arens–Dugundji 1950). Do not say cleaning is stronger than irreducibility — it is weaker (three-point example in the memo).
- **The refinement order**: AKM write $\mathfrak A\prec\mathfrak B$ for "$\mathfrak B$ refines $\mathfrak A$" (their Property 3 pins it); Sun writes $\mathcal U\succ\mathcal V$; Walkden writes $\alpha\le\beta$. The draft uses the *topology* order (refinement decreases, common refinement = meet), and the lattice meanings then fix $\vee$/$\wedge$ — DPT's point. AKM's footnote cites Kelley for terminology, but Kelley gives no symbol for cover refinement and his covering-uniformity condition is closed under *coarsening*; AKM's direction matches the **partition** convention, which is Rokhlin's ($\xi<\xi'$ iff $\xi'$ is a subpartition of $\xi$, 1967 §1.3) and older.
- **No paper found** defining $h_{\mathrm{diam}}$/$h_{\max}$ in the cover-chain sense, or stating $h_{\max}=\log\lambda_{\mathrm{Lip}}$; the closest are Kifer's pointwise characteristic exponents (metric-bound) and the adapted-metric theorems (Mather, Gourmelon — cited in `reports/diameter-growth-2cat.md` §7).

## 4. What failed, and should not be retried

- Claiming the draft's $h_{\min}$ *is* Sun's $h_L$: wrong, they sit in different columns (§2.1). Say "comparison in a fixed metric for an expansive cover", never "identity".
- An early version said the draft's $h_{\mathrm{diam}}$ line "was cut" — it was **not** cut by accident: DPT says the material was formerly in the paper and was cut when the planned scope changed. `clippings.tex` is a deliberate out-take.
- Do **not** assert a comparison between $\lambda_{\exp}$ and Schreiber's $E$ beyond $\lambda_{\exp}\ge e^E$; the sources do not compare the two classes (all compatible vs Riemannian metrics).

## 5. Open questions on this thread

1. $\lambda_{\exp}\le e^{E}$ — equivalent to $h_{\min}\le E(f,K)$ by the draft's theorem. Nothing found either way.
2. $h_{\max}=\log\lambda_{\mathrm{Lip}}$: the draft has $\le$; the converse is an adapted-metric statement, false without domination. A counterexample would live in the gap between $\sup_x\Lambda^+(x)$ and $R(f)=\lim_n\frac1n\log\sup_x\|Df^n_x\|$.
3. Whether to inline the Lipschitz-infimum and volume-side references into the diameter note — DPT said "keep the pointer", so **no** unless asked.
4. `refs/` holds 17 PDFs, including Kifer 1983, Morales–Thieullen–Villavicencio 2018, Milnor §13, Braga–Sieg, Hruska (paper and thesis), Schreiber 1997, Sun 2012, Richeson–Wiseman 2008. `refs/*.pdf` is gitignored.
5. Two Doklady notes (Kolmogorov 1958, Sinai 1959) are cited but were **not read** — paywalled, not on mathnet. The caveat is in the memo.

## 6. Conventions that matter

- DPT loses track of pronouns: **name the referent** when retracting or merging.
- Credit his contributions as his ($\lambda_{\exp}\ge e^E$; the refinement-order/lattice point; the irrational-slope objection on the diameter side).
- Report math in `$...$`; files go under the project root; **the PDF must be rebuilt whenever the `.md` changes** — he asked for this explicitly.
- Never edit `notes/sofic-status.md` — a different project.
- Commit locally per coherent unit; do not push.

## 7. Identifiers

Project root: `/workspace/nucleus/work/projects/p-492c2677edaa4757b2ea139e23b5909a/repo`. Draft clone: `.scratch/expansion` at commit `73fe2f12344ade2f4e93223f5b001724a54a846f` (2026-10-18); files `covers.tex`, `min-growth.tex`, `clippings.tex`, `alternatives.tex`, `comb-models.tex`, `visual-metrics.tex`, `intro.tex`, `Expansion.tex`, `preamble.tex`. Literature scratch: `.scratch/expansion-lit/`. Build scratch: `.scratch/report-build/<name>/`. Key commits: `da0ac2c` split the diameter note out; `f2c90dc`, `5ecddaf` bibliography trimmed to exactly what is cited; `6716cfa` 2Cat correction; later commits add the ball-cover harness and these handoffs.

## 8. Next steps, in order

1. Read `reports/expansion-prior-work.md` §1–§2 and §8 before touching anything; the memo is the state.
2. If DPT asks for more literature: the bibliographies are §9 themes A–J; the "what I did not find" lists (§6, §8.7) are the honest boundary of the search.
3. If he asks about the maximum end: the *computation* lives in `reports/diameter-growth-2cat.md`, the *literature* in §8 — keep them separate, that split was his instruction.
