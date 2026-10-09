"""Maximum-likelihood fit of the per-cycle AE noise sigma to the authors' single-copy data.

usage: python3 fit_sigma.py EXP_CQ.npy SCAN_DIR OUT_CSV [MEAS_SD]   (plot_scan.py reads fit_sigma_m0.03.csv)
Data: Cq of positive wells from ExtremeDiluteR.xlsx (github.com/Azuresky99/quPCR),
single-copy wells = Cq in [38.07, 39.5] (lower cut as in the paper; upper cut drops the
late tail the authors attribute to primer/probe mutations).
Model for each sigma: simulated Cq density (noise_scan.jl, AE = 0.9593, 5e5 runs),
convolved with Gaussian measurement noise (sd 0.03, the paper's high-copy Cq SD),
shifted by a free offset (absorbs quCq), truncated and renormalised to the window.
Profile likelihood in sigma; 95% bound from Delta(-2 ln L) = 3.84.
"""
import sys, os, re, glob
import numpy as np

EXP, SCAN, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
LO, HI = 38.07, 39.5
MEAS = float(sys.argv[4]) if len(sys.argv) > 4 else 0.03
x = np.load(EXP); x = x[(x >= LO) & (x <= HI)]
print(f"single-copy wells in [{LO}, {HI}]: {len(x)};  measurement sd = {MEAS}")
bw = 0.005
grid = np.arange(36.5, 42.0, bw)
g = np.arange(-6 * MEAS, 6 * MEAS + bw / 2, bw) if MEAS > 0 else np.array([0.0])
kern = np.exp(-0.5 * (g / max(MEAS, 1e-9)) ** 2); kern /= kern.sum()

def density(cq):
    n, _ = np.histogram(cq, np.append(grid, grid[-1] + bw))
    d = np.convolve(n.astype(float), kern, mode="same") + 1e-3   # floor: avoid log 0
    return d / (d.sum() * bw)

def loglik(d, shift):
    xs = x - shift
    inwin = (grid >= LO - shift) & (grid <= HI - shift)
    Z = d[inwin].sum() * bw
    return np.sum(np.log(np.interp(xs, grid, d) / Z))

rows = ["late,sigma,best_shift,neg2lnL"]
res = {}
for late in (0, 1):
    files = sorted((float(re.search(r"_sig([0-9.]+)\.bin", f).group(1)), f)
                   for f in glob.glob(os.path.join(SCAN, f"scan_ae0.959_late{late}_sig*.bin")))
    out = []
    for sig, f in files:
        d = density(np.fromfile(f))
        shifts = np.arange(-0.3, 0.15001, 0.0025)
        ll = np.array([loglik(d, s) for s in shifts])
        out.append((sig, shifts[ll.argmax()], -2 * ll.max()))
    m = min(o[2] for o in out)
    res[late] = out
    print(f"\nnoise {'in every cycle' if late else 'in cycles 1-20 only'}:")
    for sig, sh, v in out:
        print(f"  sigma={sig:.4f}  shift={sh:+.4f}  Delta(-2lnL)={v - m:6.2f}")
        rows.append(f"{late},{sig},{sh:.4f},{v:.3f}")
    sig = np.array([o[0] for o in out]); dv = np.array([o[2] for o in out]) - m
    best = sig[dv.argmin()]
    above = sig[(sig > best) & (dv > 3.84)]
    ub = np.interp(3.84, dv[sig >= best], sig[sig >= best]) if len(above) else np.nan
    print(f"  best sigma = {best:.4f};  95% upper bound ~ {ub:.4f}")
open(OUT, "w").write("\n".join(rows) + "\n")
