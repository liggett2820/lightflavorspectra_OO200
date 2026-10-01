#!/usr/bin/env python3
"""
Step 2b for O+O 200 GeV: rank all NBD scan points and print the centrality table of the best one.

Replaces getBestChi2_fromCat/*.sh for this system. Those scripts read fixed line numbers
("head -206") and fixed character columns ("cut -c 67-71") of the logs, which shift as soon as
one log line or number width changes. This reads the summary line that StNbdFitMaker::Scan
prints for every point,
  StNbdFitMaker::Scan  Fitting for (npp, k, x, d, chi2/ndf) = (2.400, 2.000, 0.130, 0.000, 123.4/ 80=1.543) ...
and attaches the 32 "mMultiplicityCut[...]" lines CalculateCentrality printed for that same point
(they come right before its summary line in the same log).

Usage (from the star_glauber directory, after the scan jobs finished):
  python3 find_best_nbd_OO200.py [LOG_Scan] [N best to list, default 15]
"""
import glob, os, re, sys

logdir = sys.argv[1] if len(sys.argv) > 1 else "LOG_Scan"
ntop   = int(sys.argv[2]) if len(sys.argv) > 2 else 15

summary = re.compile(r"StNbdFitMaker::Scan\s+Fitting for \(npp, k, x, d, chi2/ndf\) = \(\s*([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.eE+-]+)/\s*(\d+)=\s*([\d.eE+-]+|nan|inf)\)")
cutline = re.compile(r"(mMultiplicityCut\[[012]\]\.push_back\(.*)$")

points = []
files = sorted(glob.glob(os.path.join(logdir, "doScanX_OO200*.out")))
if not files:
    sys.exit(f"No doScanX_OO200*.out in {logdir}")

for fn in files:
    block = []
    with open(fn, errors="replace") as f:
        for line in f:
            m = cutline.search(line)
            if m:
                block.append(m.group(1).strip())
                continue
            s = summary.search(line)
            if s:
                npp, k, x, d, chi2, ndf, red = s.groups()
                try:
                    red = float(red)
                except ValueError:
                    red = float("inf")
                points.append(dict(npp=float(npp), k=float(k), x=float(x), d=float(d),
                                   chi2=float(chi2), ndf=int(ndf), chi2ndf=red, log=fn, cuts=block))
                block = []

if not points:
    sys.exit("Found logs but no scan summary lines: check that the jobs finished.")

points.sort(key=lambda p: p["chi2ndf"])
print(f"{len(points)} scan points from {len(files)} logs\n")
print(f"{'rank':>4} {'npp':>6} {'k':>6} {'x':>6} {'d':>6} {'chi2':>10} {'ndf':>4} {'chi2/ndf':>9}")
for i, p in enumerate(points[:ntop]):
    print(f"{i+1:>4} {p['npp']:6.3f} {p['k']:6.3f} {p['x']:6.3f} {p['d']:6.3f} {p['chi2']:10.1f} {p['ndf']:>4} {p['chi2ndf']:9.3f}")

# Warn if the best point sits on the edge of the scanned grid (scan range too small)
best = points[0]
for key in ("npp", "k", "x", "d"):
    vals = sorted({p[key] for p in points})
    if len(vals) > 1 and best[key] in (vals[0], vals[-1]):
        print(f"\nWARNING: best {key} = {best[key]:.3f} is at the edge of the scanned range "
              f"[{vals[0]:.3f}, {vals[-1]:.3f}]; extend the range and rescan.")

print(f"\nBest point log: {best['log']}")
print(f"Ratio file   : RatioChi2Files/Ratio_npp{best['npp']:.3f}_k{best['k']:.3f}_x{best['x']:.3f}_eff{best['d']:.3f}.root")
print("\n// ---- paste between the PASTE markers in StCentrality::Init_OO200GeV() ----")
print(f"// npp = {best['npp']:.3f}, k = {best['k']:.3f}, x = {best['x']:.3f}, d = {best['d']:.3f}, "
      f"chi2/ndf = {best['chi2']:.1f}/{best['ndf']} = {best['chi2ndf']:.3f}")
if len(best["cuts"]) != 32:
    print(f"// WARNING: expected 32 centrality lines for this point, found {len(best['cuts'])}")
for c in best["cuts"]:
    print("\t" + c)
