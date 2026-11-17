"""The extremal pairs attaining the sharp bi-Lipschitz constants for  sup |F(x)-F(y)| / |x-y|^alpha  on the Koch arc,
alpha = log3/log4, drawn on the bottom arc of the closed Koch snowflake.

    min : t = (5/16, 11/16)        chord 2/9 = 0.2222          rho = 2/6^alpha   = 0.4835
    max : t = (63/128, 321/512)    chord 2*sqrt(633)/243 = 0.2071   rho = 2*sqrt(211)/69^alpha = 1.0137

(The max is shown in its reduced (trimmed) form only.)
Run from the repository root:  python3 scripts/plot_koch_extremals.py
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

w = np.exp(1j*np.pi/3)
d = [1+0j, w, np.conj(w), 1+0j]
A = np.array([dk/3 for dk in d])
c = [0+0j]
for k in range(3):
    c.append(c[-1]+d[k])
B = np.array(c)/3.0

def F(x, depth=24):
    x = np.asarray(x, float); sh = x.shape; x = x.ravel(); n = x.size
    y = x.copy(); dig = np.zeros((depth, n), int)
    for i in range(depth):
        t = y*4; k = np.minimum(t.astype(int), 3); dig[i] = k; y = t-k
    z = np.zeros(n, complex)
    for i in range(depth-1, -1, -1):
        z = A[dig[i]]*z + B[dig[i]]
    return z.reshape(sh)

N = 9; g = np.arange(4**N+1)/4**N; K = F(g)
al = np.exp(1j*np.pi/3)
sideA, sideB, sideC = al*K, al+(1-al)*K, 1-K      # clockwise triangle, bumps outward

tmin = (5/16, 11/16)
tmax = (63/128, 321/512)                          # reduced form
P = lambda t: F(np.array(t))
pm, qm = 1-P(tmin), 1-P(tmax)
apex = 1-(0.5+1j*np.sqrt(3)/6)

fig = plt.figure(figsize=(11.0, 6.1))
ax  = fig.add_axes([0.055, 0.185, 0.545, 0.755])
axp = fig.add_axes([0.055, 0.045, 0.545, 0.085])
axi = fig.add_axes([0.625, 0.575, 0.365, 0.365])
axt = fig.add_axes([0.632, 0.10, 0.360, 0.375]); axt.axis("off")

x0, x1, y0, y1 = 0.29, 0.68, -0.335, 0.030

ax.plot(sideC.real, sideC.imag, color="0.55", lw=1.2, zorder=1)
for (a, b), col in [(tmin, "tab:blue"), (tmax, "crimson")]:
    m = (g >= a) & (g <= b)
    ax.plot((1-K[m]).real, (1-K[m]).imag, color=col, lw=3.0, alpha=0.30, zorder=2, solid_capstyle="round")
ax.plot([pm[0].real, pm[1].real], [pm[0].imag, pm[1].imag],
        color="tab:blue", lw=2.3, marker="o", ms=6, zorder=6)
ax.plot([qm[0].real, qm[1].real], [qm[0].imag, qm[1].imag],
        color="crimson", lw=2.3, marker="o", ms=6, zorder=6)
ax.plot([apex.real], [apex.imag], "k^", ms=6, zorder=7)
ax.text(apex.real+0.006, apex.imag-0.030, "$F(1/2)$", fontsize=9)
ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal"); ax.tick_params(labelsize=8)
ax.set_title(r"extremal pairs attaining the sharp bounds, on the bottom arc", fontsize=10)

for s in (sideA, sideB):
    axi.plot(s.real, s.imag, color="0.6", lw=0.9)
axi.plot(sideC.real, sideC.imag, color="k", lw=1.4)
axi.add_patch(Rectangle((x0, y0), x1-x0, y1-y0, fill=False, ec="crimson", lw=0.9, ls=(0, (4, 2))))
axi.set_xlim(0, 1.0); axi.set_ylim(-0.40, 0.95); axi.set_aspect("equal")
axi.set_xticks([]); axi.set_yticks([]); axi.set_title("closed snowflake", fontsize=8.5)
for sp in axi.spines.values():
    sp.set_color("0.85")

txt = [("tab:blue", r"$\min$:  $t=(5/16,\ 11/16)$"),
       ("tab:blue", r"$d_K=2/9=0.2222$, gap $3/8$"),
       ("tab:blue", r"$\rho=2/6^{\alpha}=0.4835$"),
       ("crimson", r"$\max$:  $t=(63/128,\ 321/512)$"),
       ("crimson", r"$d_K=2\sqrt{633}/243=0.2071$, gap $69/512$"),
       ("crimson", r"$\rho=2\sqrt{211}/69^{\alpha}=1.0137$"),
       ("0.35",    r"$\rho=d_K/|x-y|^{\alpha}$")]
y = 0.97
for col, s in txt:
    axt.text(0, y, s, fontsize=9, color=col, transform=axt.transAxes); y -= 0.135

axp.axhline(0, color="0.6", lw=1)
for xv, c_, dy, lab in [(tmin[0], "tab:blue", .38, r"$5/16$"), (tmin[1], "tab:blue", .38, r"$11/16$"),
                        (tmax[0], "crimson", -.52, r"$63/128$"), (tmax[1], "crimson", -.52, r"$321/512$")]:
    axp.plot([xv, xv], [0, dy*0.45], color=c_, lw=1)
    axp.text(xv, dy, lab, color=c_, fontsize=8.2, ha="center", va="center")
axp.plot([0.5], [0], "k^", ms=5); axp.text(0.5, 0.38, r"$t=1/2$", fontsize=8.2, ha="center")
axp.set_ylim(-0.68, 0.62); axp.set_xlim(0.81, 0.22); axp.axis("off")
axp.set_title(r"parameter $t$ (increasing leftwards, matching the drawing)", fontsize=8.2, pad=0)
fig.savefig("figures/koch-extremal-configurations.png", dpi=170, bbox_inches="tight", facecolor="white")
print("saved")
