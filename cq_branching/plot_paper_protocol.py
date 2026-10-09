"""Plot the authors'-protocol runs (paper_protocol.jl) with their binning and smoothing.

usage: python3 plot_paper_protocol.py DATA_DIR OUT_DIR
Histogram: 1000 equal bins over the full data range (MATLAB hist(Cq,1000)).
Smoothing: Gaussian template of width W bins, sd W/5, where W is 0.15 cycles
(Fig. 1D / S3A script) or one tenth of the 5-95% range (Fig. S3B script).
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA, OUT = sys.argv[1], sys.argv[2]
BLUE, RED, INK, MUTED, GOLD = "#3F6FD1", "#D9442B", "#1F2430", "#6B7280", "#E0A526"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

def matlab_smooth(cq, mode):
    counts, edges = np.histogram(cq, bins=1000)
    centers = 0.5 * (edges[1:] + edges[:-1]); dx = centers[1] - centers[0]
    if mode == "fixed":
        w = round(0.15 / dx)
    else:
        cp = np.cumsum(counts) / counts.sum()
        i5, i95 = np.argmax(cp >= 0.05), np.argmax(cp >= 0.95)
        w = round((centers[i95] - centers[i5]) / 10 / dx)
    w = (w // 2) * 2 + 1; h = w // 2
    k = np.exp(-0.5 * ((np.arange(1, w + 1) - (h + 1)) / (w / 5)) ** 2) / (np.sqrt(2 * np.pi) * w / 5)
    sm = counts.astype(float).copy()
    for i in range(h, len(counts) - h):
        sm[i] = np.sum(counts[i - h:i + h + 1] * k)
    sl = slice(h, len(counts) - h)
    return centers, counts, centers[sl], sm[sl], dx

runs = [("paper_fixed.bin", "fixed", "AE = 0.959 (authors' protocol)", 38.21),
        ("paper_noise.bin", "noise", "AE = 0.959 ± 0.1 (95% half-width)\nnoise in cycles 1-20 only", 38.31),
        ("paper_noise_all.bin", "noise", "AE = 0.959 ± 0.1 (95% half-width)\nnoise in every cycle", None)]
fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.8))
for ax, (fn, mode, lab, paper), letter in zip(axes, runs, "abc"):
    cq = np.fromfile(os.path.join(DATA, fn), dtype="<f8")
    c, n, cs, sm, dx = matlab_smooth(cq, mode)
    ax.bar(c, n, width=dx, color=BLUE, linewidth=0, label="Cq")
    ax.plot(cs, sm, color=RED, lw=1.2, label="smooth Cq")
    cqp = cs[np.argmax(sm)]
    ax.axvline(38.34, color=GOLD, lw=1.2)
    txt = f"{lab}\nCqP = {cqp:.2f}" + (f"  (paper {paper:.2f})" if paper else "") + f"\nstd = {cq.std():.3f}"
    ax.text(0.97, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=9, color=INK, linespacing=1.4)
    ax.set_xlim(38, 40); ax.set_ylim(0, max(n[(c > 38) & (c < 40)].max(), sm.max()) * 1.45)
    ax.set_xlabel("Cq"); ax.set_ylabel("Frequency")
    ax.grid(axis="y", color="#E5E7EB", lw=0.6); ax.set_axisbelow(True)
    ax.text(-0.12, 1.04, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", color=INK)
    print(f"{fn:22s} CqP={cqp:.3f} mean={cq.mean():.3f} std={cq.std():.3f}")
axes[0].legend(loc="center right", frameon=False, fontsize=9)
fig.suptitle("Authors' protocol: 19 stochastic cycles, then deterministic growth to (1.959)^38.34; 1e6 runs; yellow line = quCq 38.34",
             fontsize=10, color=MUTED, y=1.0)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_paper_protocol.png"), dpi=160, bbox_inches="tight")
