"""Peak survival vs per-cycle AE noise (output of noise_scan.jl).

usage: python3 analyze_scan.py SCAN_DIR OUT_DIR AE [AE ...]
Histogram: 0.005-cycle bins; smoothing: Gaussian, 0.15-cycle window with sd = window/5
(the authors' experimental smoothing). Side peaks are located on the sigma = 0 curve.
Visibility of side peak j: V_j = (max_j - min_j) / max_j, where max_j is the local maximum
nearest the sigma = 0 position (within 0.06 cycles) and min_j the minimum between it
and the previous peak; V_j = 0 if no local maximum exists there.
Writes visibility_<AE>.csv and prints the table.
"""
import sys, os, re, glob
import numpy as np

SCAN, OUT, AES = sys.argv[1], sys.argv[2], [float(a) for a in sys.argv[3:]]
BW = 0.005
W = round(0.15 / BW); W = (W // 2) * 2 + 1; H = W // 2
K = np.exp(-0.5 * ((np.arange(W) - H) / (W / 5)) ** 2); K /= K.sum()

def smooth(cq):
    lo, hi = np.floor(cq.min()) - 0.5, np.ceil(np.percentile(cq, 99.9)) + 0.5
    e = np.arange(lo, hi + BW, BW)
    n, _ = np.histogram(cq, e)
    return 0.5 * (e[1:] + e[:-1]), np.convolve(n, K, mode="same")

def maxima(c, s):
    i = np.where((s[1:-1] > s[:-2]) & (s[1:-1] >= s[2:]) & (s[1:-1] > 0.003 * s.max()))[0] + 1
    return c[i], s[i]

def load(ae, late):
    rows = []
    for f in glob.glob(os.path.join(SCAN, f"scan_ae{ae}_late{late}_sig*.bin")):
        rows.append((float(re.search(r"_sig([0-9.]+)\.bin", f).group(1)), f))
    return sorted(rows)

for ae in AES:
    ref_c, ref_s = smooth(np.fromfile(load(ae, 0)[0][1]))
    mc, ms = maxima(ref_c, ref_s)
    x1 = ref_c[np.argmax(ref_s)]
    ref_peaks = [x for x in mc if x <= x1 + 1.6]          # peak 1 and its side peaks
    print(f"\nAE = {ae}: sigma = 0 maxima at {np.round(mc, 3).tolist()}  (offsets {np.round(np.array(ref_peaks) - x1, 3).tolist()})")
    lines = ["late,sigma," + ",".join(f"V@{x - x1:+.2f}" for x in ref_peaks[1:]) + ",cqp,std,n_maxima"]
    for late in (0, 1):
        print(f"  noise {'in every cycle' if late else 'in cycles 1-20 only'}:")
        for sig, f in load(ae, late):
            cq = np.fromfile(f)
            c, s = smooth(cq)
            mx, my = maxima(c, s)
            p1 = c[np.argmax(s)]
            shift = p1 - x1
            V = []
            prev = p1
            for x in ref_peaks[1:]:
                target = x + shift
                j = np.where(np.abs(mx - target) < 0.06)[0]
                if len(j) == 0:
                    V.append(0.0); continue
                xm, ym = mx[j[0]], my[j[0]]
                seg = (c > prev) & (c < xm)
                ymin = s[seg].min() if seg.any() else ym
                V.append(max(0.0, (ym - ymin) / ym)); prev = xm
            nmax = int(np.sum((mx > p1 - 0.5) & (mx < p1 + 1.6)))
            lines.append(f"{late},{sig},{','.join(f'{v:.3f}' for v in V)},{p1:.3f},{cq.std():.3f},{nmax}")
            print(f"    sigma={sig:.4f}  CqP={p1:.3f}  std={cq.std():.3f}  maxima={nmax}  V={[round(v, 3) for v in V]}")
    open(os.path.join(OUT, f"visibility_{ae}.csv"), "w").write("\n".join(lines) + "\n")
