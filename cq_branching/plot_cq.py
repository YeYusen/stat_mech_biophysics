"""Static figures from the Julia Cq simulation (sim.jl).

usage: python3 plot_cq.py DATA_DIR OUT_DIR
Reads cq_<mu>_<sigma>[_extend].bin (Float64 arrays) and writes
  fig_cq_0.959.png  -- reproduction of the paper's Fig. S3a/b (+ p>1 "extend" variant)
  fig_cq_0.8.png    -- same for AE = 0.8
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

DATA, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

BLUE, RED, INK, MUTED = "#3F6FD1", "#D9442B", "#1F2430", "#6B7280"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK, "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False,
                     "axes.spines.right": False})

def load(name):
    return np.fromfile(os.path.join(DATA, name), dtype="<f8")

def delta_k(mu, k):
    """Cq delay from a single failed duplication at cycle k (one of 2^(k-1) particles)."""
    return np.log2(2.0**k / (2.0**k - 1.0)) / np.log2(1.0 + mu)

def panel(ax, cq, mu, label, xlim, bw=0.01, annotate_peaks=True, letter=""):
    bins = np.arange(xlim[0], xlim[1] + bw, bw)
    h, e = np.histogram(cq, bins=bins)
    c = 0.5 * (e[1:] + e[:-1])
    ax.bar(c, h, width=bw, color=BLUE, linewidth=0, label="Cq")
    sm = gaussian_filter1d(h.astype(float), sigma=2.0)
    ax.plot(c, sm, color=RED, lw=1.2, label="smooth Cq")
    mode = c[np.argmax(sm)]
    mean = cq.mean()
    ax.axvline(mean, color="#E0A526", lw=1.2)
    if annotate_peaks:
        ax.annotate("1", xy=(mode, sm.max()),
                    xytext=(mode - 0.09 * (xlim[1] - xlim[0]), sm.max() * 1.02),
                    color=INK, fontsize=9, ha="center", va="bottom",
                    arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
        # side peaks: one failed duplication at cycle k = 3, 2, 1  ->  labels 2, 3, 4
        for lab, k in zip((2, 3, 4), (3, 2, 1)):
            x = mode + delta_k(mu, k)
            i = np.argmin(np.abs(c - x))
            y = sm[max(i - 5, 0): i + 6].max()
            ax.annotate(f"{lab}", xy=(x, y), xytext=(x + 0.04 * (xlim[1] - xlim[0]), y + 0.10 * sm.max()),
                        color=INK, fontsize=9, ha="center",
                        arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
    info = (f"{label}\nCqP (mode) = {mode:.2f}\nmean = {mean:.2f} (yellow line)\nstd = {cq.std():.3f}")
    ax.text(0.97, 0.97, info, transform=ax.transAxes, ha="right", va="top", color=INK, fontsize=9,
            linespacing=1.4)
    ax.set_xlim(*xlim)
    ax.set_ylim(0, sm.max() * 1.18)
    ax.set_xlabel("Cq")
    ax.set_ylabel("Frequency")
    ax.grid(axis="y", color="#E5E7EB", lw=0.6)
    ax.set_axisbelow(True)
    if letter:
        ax.text(-0.12, 1.04, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", color=INK)
    return mode

def figure(mu, xlim, fname, title):
    cq0 = load(f"cq_{mu}_0.0.bin")
    cq1 = load(f"cq_{mu}_0.1.bin")
    cq2 = load(f"cq_{mu}_0.1_extend.bin")
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.8))
    m0 = panel(axes[0], cq0, mu, f"AE = {mu}", xlim, letter="a")
    panel(axes[1], cq1, mu, f"AE = {mu} ± 0.1, p > 1 clipped to 1", xlim, annotate_peaks=False, letter="b")
    panel(axes[2], cq2, mu, f"AE = {mu} ± 0.1, p > 1 gives extra copies", xlim, annotate_peaks=False, letter="c")
    axes[0].legend(loc="upper right", frameon=False, fontsize=9, bbox_to_anchor=(1.0, 0.60))
    fig.suptitle(title, fontsize=10, color=MUTED, y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=160, bbox_inches="tight")
    plt.close(fig)
    # print summary
    for name, cq in (("sigma=0", cq0), ("sigma=0.1 clip", cq1), ("sigma=0.1 extend", cq2)):
        print(f"mu={mu} {name:18s} n={len(cq)} mean={cq.mean():.3f} std={cq.std():.3f}")
    print(f"mu={mu} predicted side-peak offsets from main peak: "
          + ", ".join(f"k={k}: +{delta_k(mu,k):.3f}" for k in (3, 2, 1)))
    return m0

figure(0.959, (38.7, 40.7), "fig_cq_0.959.png",
       "Branching-process PCR, N0 = 1, threshold 2^38, 2e5 runs each; log-interpolated Cq")
figure(0.8, (43.3, 47.3), "fig_cq_0.8.png",
       "Branching-process PCR, N0 = 1, threshold 2^38, 2e5 runs each; log-interpolated Cq")
