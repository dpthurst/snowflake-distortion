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

`ladder-m.png` — best known configurations at m=100,200,400,800,1600.
`ladder-D.png` — D_m against m, with the explicit bounds as reference lines.
Figure 3 of notes/improving-koch-extremals.pdf.
