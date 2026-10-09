"""Wells needed to bound the per-cycle AE noise, from the simulated scans.

usage: python3 wells_needed.py SCAN_DIR OUT_CSV
If the truth is sigma = 0, the expected likelihood-ratio statistic against sigma with N
single-copy wells is 2 N KL(p_0 || p_sigma) (shift profiled out). The 95% bound is reached
when this equals 3.84, so N_95(sigma) = 1.92 / min_shift KL. Densities: 0.01-cycle bins,
Gaussian measurement noise sd 0.03, central 99% of the sigma = 0 distribution as window.
"""
import sys, os, re, glob
import numpy as np

SCAN, OUT = sys.argv[1], sys.argv[2]
MEAS, BW = 0.03, 0.01

def load(ae, late):
    rows = []
    for f in glob.glob(os.path.join(SCAN, f"scan_ae{ae}_late{late}_sig*.bin")):
        rows.append((float(re.search(r"_sig([0-9.]+)\.bin", f).group(1)), f))
    return sorted(rows)

g = np.arange(-5 * MEAS, 5 * MEAS + BW / 2, BW); kern = np.exp(-0.5 * (g / MEAS) ** 2); kern /= kern.sum()
lines = ["ae,late,sigma,KL,N95"]
for ae in (0.959, 0.8):
    for late in (0, 1):
        rows = load(ae, late)
        c0 = np.fromfile(rows[0][1])
        lo, hi = np.percentile(c0, [0.5, 99.5])
        edges = np.arange(lo - 1.0, hi + 1.0, BW)
        def dens(cq):
            n, _ = np.histogram(cq, edges)
            d = np.convolve(n.astype(float), kern, mode="same") + 0.5
            return d
        d0 = dens(c0); cen = 0.5 * (edges[1:] + edges[:-1])
        win = (cen >= lo) & (cen <= hi)
        p0 = d0[win] / d0[win].sum()
        print(f"\nAE={ae}  noise {'every cycle' if late else 'cycles 1-20'}  window [{lo:.2f}, {hi:.2f}]")
        for sig, f in rows[1:]:
            d = dens(np.fromfile(f))
            best = np.inf
            for k in range(-40, 41):                       # shift by k bins (+-0.4 cycles)
                ds = np.roll(d, k)[win]; q = ds / ds.sum()
                best = min(best, np.sum(p0 * np.log(p0 / q)))
            n95 = 1.92 / best
            lines.append(f"{ae},{late},{sig},{best:.6g},{n95:.4g}")
            if abs(sig * 400 - round(sig * 400)) < 1e-9 and round(sig * 400) % 2 == 0:
                print(f"  sigma={sig:.3f}  KL={best:.2e}  N95={n95:9.0f}")
open(OUT, "w").write("\n".join(lines) + "\n")
