# figures

`koch-extremal-configurations.png` — the two extremal pairs for
`sup |F(x)-F(y)| / |x-y|^alpha` on the Koch arc (alpha = log3/log4), drawn on the
bottom arc of the closed Koch snowflake. They are no longer conjectural: both
constants are proved and attained, see notes/koch-extremals.tex. The figure was
drawn while they were still conjectures.
min: t=(5/16, 11/16), chord 2/9 = 0.2222
max: t=(63/256, 321/1024), chord 2*sqrt(211)/243 = 0.1196
(trimmed representative of the same max: t=(63/128, 321/512), chord 0.2071)

Reproduce with:  python3 scripts/plot_koch_extremals.py

`distortion-heatmaps.png` — r(x,y) = |F(x)-F(y)|/|x-y|^alpha on the parameter square for the
Koch arc, the two-piece and three-piece reversible arcs (each on its own parameter grid) and the
m=800 ladder configuration; each panel normalised to its own range, blue = min, white = geometric
mean, red = max, with the Koch extremal pairs marked. Figure 3 of the note.

`small-m-optima.png` — the optimal configurations at m=3,4,5, endpoint-gauged with the nodes
numbered: a single corner (D=1 exactly), an M (1.0628) and a W (1.0931). Figure 4 of the note.

`block-similarity.png` — the best matches between different parts of the m=1600 configuration, one
per bin of size ratio. Top: block solid, the matching arc dashed after normalising and turning it
back (residuals 0.10-0.17 at size ratios 0.51-1.78). Bottom: where the two pieces sit on the whole
curve, highlighted. Figure 5 of the note.

`block-comparison.png` — deviation of the restricted distortion of every contiguous block of the
m=800 configuration from the whole, against block length, with the exact level a self-similar arc
attains at its substitution scales. Figure 4 of the note.

`explicit-curves.png` — Koch, the two-piece reversible arc and the three-piece
reversible arc (the two explicit arcs below the Koch constant). Reproduce with:
python3 scripts/plot_explicit_and_ladder.py

`ladder-m.png` — best known S-branch configurations at m=100,200,400,800,1600 (dotted:
the chord, which the S crosses).
`ladder-D.png` — D_m against m between the certified bounds (purple: D_4, red: Koch), with
the bar to the other branch at each m, and that gap as a percentage. Figure 5 of the note.

`uneven-spacing.png` — finite-m optima with the parameter gaps bounded: best D_m against m for
gap ratios 1..32, the m=16 shapes overlaid, and the optimal spacings at m=32.
`uneven-random-shapes.png` — random fixed spacings at m=16,32,64 against the m=800 configuration.
`uneven-random-diagnostics.png` — D against m for three random draws against the uniform ladder,
the three m=64 draws overlaid, and their distance to the m=800 curve.
`stress-measures.png` — the max-entropy stress at m=200,400,800: the pair square, the marginal of
the min measure, the min measure in lattice units and the max measure against separation.
Preliminary, not in either note; see §13 of results/uneven-spacing-preliminary.md.

`uneven-ladder.png` — the bounded-κ ladder: D against m for κ = 1,2,4,8,16,32 (κ=4 extended to
m=400) with the polished uniform ladder, and the increment per doubling against m.
All of these are preliminary and deliberately not in either note; write-up in
results/uneven-spacing-preliminary.md.

`arc-m.png` — the one-sided (arc) branch at m=40..1600: the gauged curve touches its chord
but never crosses it. Figure 4 of the note.

`construction-stages.png` — the two explicit constructions at levels 1,2,3 and at the
level where the picture settles: each segment is replaced by N scaled/rotated copies
of the whole, some reflected and some traversed backwards. Titles give D measured on
the drawn polyline (1.00000, 1.89800, 1.89800, 2.04147 for two pieces;
1.07457, 1.98927, 2.01124, 2.01553 for three). Figure 2 of notes/improving-koch-extremals.pdf.
Reproduce with: python3 scripts/plot_construction_stages.py
