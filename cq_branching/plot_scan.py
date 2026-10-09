"""Figures for the noise scan in the authors' protocol and the fit to their data.

usage: python3 plot_scan.py SCAN_DIR DATA_DIR EXP_CQ.npy OUT_DIR
  fig_scan_<AE>.png       histograms at selected sigma, rows: noise in cycles 1-20 / every cycle
  fig_visibility.png      side-peak visibility vs sigma, AE = 0.959 (visibility_0.959.csv)
  fig_wells_needed.png    wells needed for a 95% bound (wells_needed.csv)
  fig_fit_experiment.png  experimental single-copy Cq vs model, and the profile likelihood
"""
import sys, os, re, glob, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCAN, DATA, EXP, OUT = sys.argv[1:5]
BLUE, RED, INK, MUTED, GOLD = "#3F6FD1", "#D9442B", "#1F2430", "#6B7280", "#E0A526"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
BW = 0.005
W = round(0.15 / BW); W = (W // 2) * 2 + 1; H = W // 2
K = np.exp(-0.5 * ((np.arange(W) - H) / (W / 5)) ** 2); K /= K.sum()

def f_for(ae, late, sig):
    return os.path.join(SCAN, f"scan_ae{ae}_late{late}_sig{round(sig, 4)}.bin")

def hist_smooth(cq, lo, hi, bw=BW):
    e = np.arange(lo, hi + bw, bw)
    n, _ = np.histogram(cq, e)
    return 0.5 * (e[1:] + e[:-1]), n, np.convolve(n, K, mode="same")

PICK = [0.0, 0.005, 0.01, 0.02, 0.035, 0.05]
for ae, xlim in ((0.959, (38.0, 40.0)), (0.8, (42.8, 46.6))):
    if not os.path.exists(f_for(ae, 1, 0.05)):
        print("missing scan for", ae); continue
    fig, axes = plt.subplots(2, len(PICK), figsize=(3.0 * len(PICK), 5.4), sharey="row")
    for row, late in enumerate((0, 1)):
        ref = None
        for col, sig in enumerate(PICK):
            cq = np.fromfile(f_for(ae, late, sig))
            c, n, s = hist_smooth(cq, xlim[0] - 0.3, xlim[1] + 0.3)
            ax = axes[row, col]
            ax.bar(c, n, width=BW, color=BLUE, linewidth=0)
            ax.plot(c, s, color=RED, lw=1.0)
            if ref is None: ref = s.max() * 1.1
            ax.set_xlim(*xlim); ax.set_ylim(0, ref)
            ax.set_title(f"σ = {sig:g}", fontsize=10, color=INK)
            ax.grid(axis="y", color="#E5E7EB", lw=0.6); ax.set_axisbelow(True)
            if col == 0:
                ax.set_ylabel(("noise in cycles 1-20\n" if late == 0 else "noise in every cycle\n") + "Frequency")
            if row == 1: ax.set_xlabel("Cq")
    fig.suptitle(f"AE = {ae}, authors' protocol, threshold 1.959^38.34, 5e5 runs per panel; "
                 f"per-cycle AE ~ N({ae}, σ²)", fontsize=10, color=MUTED)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_scan_{ae}.png"), dpi=140)
    plt.close(fig)

# visibility (AE = 0.959; at AE = 0.8 there are no separate side peaks even at sigma = 0)
fig, ax = plt.subplots(1, 1, figsize=(6.5, 3.6))
rows = list(csv.reader(open(os.path.join(DATA, "visibility_0.959.csv"))))
hdr, rows = rows[0], rows[1:]
vcols = [i for i, h in enumerate(hdr) if h.startswith("V@")]
names = {"+0.39": "peak 3 (+0.39)", "+1.00": "peak 4 (+1.00)", "+1.39": "peak 3 copy (+1.39)"}
colors = [RED, BLUE, GOLD]
for j, i in enumerate(vcols):
    for late, ls in ((0, "-"), (1, "--")):
        r = [x for x in rows if x[0] == str(late)]
        ax.plot([float(x[1]) for x in r], [float(x[i]) for x in r], ls, color=colors[j % 3], lw=1.6,
                marker="o", ms=3, label=f"{names.get(hdr[i][2:], hdr[i][2:])}, {'cycles 1-20' if late == 0 else 'every cycle'}")
ax.set_xlabel("per-cycle σ of AE"); ax.set_ylabel("visibility (max − dip) / max")
ax.set_title("AE = 0.959, authors' protocol", fontsize=10)
ax.grid(color="#E5E7EB", lw=0.6); ax.legend(frameon=False, fontsize=7.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_visibility.png"), dpi=150); plt.close(fig)

# wells needed for a 95% bound, from wells_needed.py
fn = os.path.join(DATA, "wells_needed.csv")
if os.path.exists(fn):
    rows = list(csv.reader(open(fn)))[1:]
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 3.6))
    for ae, col in ((0.959, BLUE), (0.8, RED)):
        for late, ls in ((0, "-"), (1, "--")):
            r = [x for x in rows if float(x[0]) == ae and x[1] == str(late)]
            ax.plot([float(x[2]) for x in r], [float(x[4]) for x in r], ls, color=col, lw=1.6, marker="o", ms=3,
                    label=f"AE = {ae}, {'cycles 1-20' if late == 0 else 'every cycle'}")
    ax.axhline(297, color=MUTED, lw=1, ls=":"); ax.text(0.0255, 330, "297 wells (experiment)", color=MUTED, fontsize=8)
    ax.axhspan(1e4, 1e6, color="#F3F4F6", zorder=0); ax.text(0.0255, 1.3e4, "Monte Carlo noise floor", color=MUTED, fontsize=8)
    ax.set_yscale("log"); ax.set_ylim(10, 1e5)
    ax.set_xlabel("per-cycle σ of AE to exclude"); ax.set_ylabel("single-copy wells for 95% exclusion")
    ax.grid(color="#E5E7EB", lw=0.6, which="both"); ax.legend(frameon=False, fontsize=7.5)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_wells_needed.png"), dpi=150); plt.close(fig)

# fit to experiment
x = np.load(EXP); xs = x[(x >= 38.07) & (x <= 39.5)]
fit = list(csv.reader(open(os.path.join(DATA, "fit_sigma_m0.03.csv"))))[1:]
fig, axes = plt.subplots(1, 2, figsize=(12, 3.8), gridspec_kw={"width_ratios": [1.5, 1]})
ax = axes[0]
bw = 0.05
e = np.arange(38.0, 39.55, bw)
ax.hist(xs, bins=e, color=BLUE, alpha=0.85, label=f"experiment, {len(xs)} single-copy wells")
for sig, col, ls in ((0.0, INK, "-"), (0.015, RED, "-"), (0.035, GOLD, "-")):
    late = 1
    r = [float(rr[2]) for rr in fit if rr[0] == str(late) and abs(float(rr[1]) - sig) < 1e-9]
    shift = r[0] if r else 0.0
    cq = np.fromfile(f_for(0.959, late, sig)) + shift
    cq = cq + np.random.default_rng(1).normal(0, 0.03, len(cq))
    cq = cq[(cq >= 38.07) & (cq <= 39.5)]
    h, _ = np.histogram(cq, e)
    ax.step(0.5 * (e[1:] + e[:-1]), h / h.sum() * len(xs), where="mid", color=col, lw=1.4, ls=ls,
            label=f"model σ = {sig:g} (every cycle)")
ax.set_xlabel("Cq"); ax.set_ylabel("wells per 0.05 cycle"); ax.legend(frameon=False, fontsize=8)
ax.set_xlim(38.0, 39.5)
ax = axes[1]
for late, col in ((0, BLUE), (1, RED)):
    r = [rr for rr in fit if rr[0] == str(late)]
    sig = np.array([float(rr[1]) for rr in r]); v = np.array([float(rr[3]) for rr in r])
    ax.plot(sig, v - v.min(), "-o", ms=3, color=col, label="noise in cycles 1-20" if late == 0 else "noise in every cycle")
ax.axhline(3.84, color=MUTED, lw=1, ls="--"); ax.text(0.036, 2.9, "95% (Δ = 3.84)", color=MUTED, fontsize=8)
ax.set_ylim(0, 20); ax.set_xlabel("per-cycle σ of AE"); ax.set_ylabel("Δ(−2 ln L)")
ax.legend(frameon=False, fontsize=8); ax.grid(color="#E5E7EB", lw=0.6)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_fit_experiment.png"), dpi=150); plt.close(fig)
print("done")
