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

`explicit-curves.png` — Koch, the two-piece reversible arc and the three-piece
reversible arc (the two explicit arcs below the Koch constant). Reproduce with:
python3 scripts/plot_explicit_and_ladder.py

`ladder-m.png` — best known S-branch configurations at m=100,200,400,800,1600 (dotted:
the chord, which the S crosses).
`ladder-D.png` — D_m against m between the certified bounds (purple: D_4, red: Koch), with
the bar to the other branch at each m, and that gap as a percentage. Figure 5 of the note.

`arc-m.png` — the one-sided (arc) branch at m=40..1600: the gauged curve touches its chord
but never crosses it. Figure 4 of the note.

`construction-stages.png` — the two explicit constructions at levels 1,2,3 and at the
level where the picture settles: each segment is replaced by N scaled/rotated copies
of the whole, some reflected and some traversed backwards. Titles give D measured on
the drawn polyline (1.00000, 1.89800, 1.89800, 2.04147 for two pieces;
1.07457, 1.98927, 2.01124, 2.01553 for three). Figure 2 of notes/improving-koch-extremals.pdf.
Reproduce with: python3 scripts/plot_construction_stages.py
