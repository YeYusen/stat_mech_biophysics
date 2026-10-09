"""Animation of the Cq histogram as the per-cycle AE noise sigma sweeps 0 -> 0.1.

usage: python3 animate_cq.py SWEEP_DIR_0.959 SWEEP_DIR_0.8 OUT_DIR
Reads the cq_mu<mu>_sig<sigma>.bin files written by sweep.jl, renders one PNG
per sigma (two panels: AE = 0.959 and AE = 0.8), then calls ffmpeg to make
cq_sigma_sweep.mp4 and cq_sigma_sweep.gif.
"""
import sys, os, re, glob, subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

D1, D2, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
FR = os.path.join(OUT, "frames"); os.makedirs(FR, exist_ok=True)
BLUE, RED, INK, MUTED = "#3F6FD1", "#D9442B", "#1F2430", "#6B7280"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK, "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False,
                     "axes.spines.right": False})

def files(d):
    out = []
    for f in sorted(glob.glob(os.path.join(d, "cq_mu*_sig*.bin"))):
        s = float(re.search(r"_sig([0-9.]+)\.bin", f).group(1))
        out.append((s, f))
    return sorted(out)

F1, F2 = files(D1), files(D2)
assert len(F1) == len(F2) and all(abs(a[0] - b[0]) < 1e-9 for a, b in zip(F1, F2)), "sigma grids differ"

SPEC = [  # (mu, xlim, bin width)
    (0.959, (38.7, 41.0), 0.01),
    (0.8,   (43.3, 47.5), 0.02),
]
# fixed y-limits from the sigma = 0 frame (tallest peaks) so the eye sees the broadening
ymax = []
for (mu, xlim, bw), (s, f) in zip(SPEC, (F1[0], F2[0])):
    cq = np.fromfile(f, dtype="<f8")
    h, _ = np.histogram(cq, bins=np.arange(xlim[0], xlim[1] + bw, bw))
    ymax.append(gaussian_filter1d(h.astype(float), 2.0).max() * 1.12)

def delta_k(mu, k):
    return np.log2(2.0**k / (2.0**k - 1.0)) / np.log2(1.0 + mu)

for i, ((s, f1), (_, f2)) in enumerate(zip(F1, F2)):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, (mu, xlim, bw), f, ym in zip(axes, SPEC, (f1, f2), ymax):
        cq = np.fromfile(f, dtype="<f8")
        bins = np.arange(xlim[0], xlim[1] + bw, bw)
        h, e = np.histogram(cq, bins=bins)
        c = 0.5 * (e[1:] + e[:-1])
        ax.bar(c, h, width=bw, color=BLUE, linewidth=0)
        sm = gaussian_filter1d(h.astype(float), 2.0)
        ax.plot(c, sm, color=RED, lw=1.2)
        ax.set_xlim(*xlim); ax.set_ylim(0, ym)
        ax.set_xlabel("Cq"); ax.set_ylabel("Frequency")
        ax.grid(axis="y", color="#E5E7EB", lw=0.6); ax.set_axisbelow(True)
        ax.text(0.97, 0.92, f"AE = {mu} ± {s:.4f}".rstrip("0").rstrip("."),
                transform=ax.transAxes, ha="right", va="top", color=INK, fontsize=11)
        ax.text(0.97, 0.80, f"mean {cq.mean():.2f}, std {cq.std():.3f}",
                transform=ax.transAxes, ha="right", va="top", color=MUTED, fontsize=9)
    fig.suptitle(f"Cq distribution vs. per-cycle AE noise   (sigma = {s:.4f};  "
                 f"1e5 runs per frame, p > 1 clipped to 1)", fontsize=10, color=MUTED)
    fig.tight_layout()
    fig.savefig(os.path.join(FR, f"frame_{i:03d}.png"), dpi=110)
    plt.close(fig)
print("frames:", len(F1))

# static filmstrip: six sigma values, both mu, for the README
PICK = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]
fig, axes = plt.subplots(2, len(PICK), figsize=(3.0 * len(PICK), 5.6))
for row, ((mu, xlim, bw), FF, ym) in enumerate(zip(SPEC, (F1, F2), ymax)):
    for col, sp in enumerate(PICK):
        s, f = min(FF, key=lambda t: abs(t[0] - sp))
        cq = np.fromfile(f, dtype="<f8")
        ax = axes[row, col]
        bins = np.arange(xlim[0], xlim[1] + bw, bw)
        h, e = np.histogram(cq, bins=bins)
        c = 0.5 * (e[1:] + e[:-1])
        ax.bar(c, h, width=bw, color=BLUE, linewidth=0)
        ax.plot(c, gaussian_filter1d(h.astype(float), 2.0), color=RED, lw=1.0)
        ax.set_xlim(*xlim); ax.set_ylim(0, ym)
        ax.set_title(f"AE = {mu} ± {s:g}", fontsize=10, color=INK)
        ax.text(0.97, 0.95, f"std {cq.std():.2f}", transform=ax.transAxes, ha="right", va="top",
                color=MUTED, fontsize=8)
        ax.grid(axis="y", color="#E5E7EB", lw=0.6); ax.set_axisbelow(True)
        if col: ax.set_yticklabels([])
        else:   ax.set_ylabel("Frequency")
        if row: ax.set_xlabel("Cq")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "cq_sigma_filmstrip.png"), dpi=130)
plt.close(fig)

mp4 = os.path.join(OUT, "cq_sigma_sweep.mp4")
gif = os.path.join(OUT, "cq_sigma_sweep.gif")
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "4",
                "-i", os.path.join(FR, "frame_%03d.png"),
                "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", mp4], check=True)
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "4",
                "-i", os.path.join(FR, "frame_%03d.png"),
                "-vf", "fps=4,scale=900:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=64[p];[s1][p]paletteuse=dither=bayer:bayer_scale=3",
                "-loop", "0", gif], check=True)
print("wrote", mp4, gif)
