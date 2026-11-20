#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qe_eos_fit.py -- fit an equation of state (Birch-Murnaghan 3rd order / Murnaghan) to E(V) or E(a).
                 由一组 (V,E) 或 (a,E) 数据拟合三阶 Birch-Murnaghan（或 Murnaghan）状态方程，
                 输出 V0 / a0、B0（GPa）、B0'，并可直接从一组 pw.x 输出中提取体积和能量。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Input
-----
  Data file: two columns, '#' comments allowed.
     --x volume  : V in A^3 (default) or bohr^3 (--xunit bohr)
     --x alat    : lattice constant a in A (or bohr), V = f * a^3 with --vfactor f
                   (fcc primitive cell f = 1/4, bcc 1/2, simple cubic 1, diamond/zincblende 1/4)
     energies in Ry (default, pw.x) or eV (--eunit eV)
  or --from-out *.out : read "unit-cell volume = ... (a.u.)^3" and the last
     "!    total energy = ... Ry" of each pw.x output (one scf per volume).
Equations (E0, V0, B0, B0p):
  BM3:  E = E0 + 9 V0 B0/16 { [(V0/V)^(2/3) - 1]^3 B0p + [(V0/V)^(2/3) - 1]^2 [6 - 4 (V0/V)^(2/3)] }
  Murnaghan: E = E0 + B0 V/B0p [ (V0/V)^B0p/(B0p - 1) + 1 ] - B0 V0/(B0p - 1)
Units: B0 reported in GPa (1 eV/A^3 = 160.21766 GPa). Energies per cell as given.
Check: points should bracket the minimum (about +-5 % in a / +-15 % in V); the script warns
if the minimum lies outside the data range.  QE's own ev.x does the same fits interactively.

Usage examples
--------------
    python qe_eos_fit.py ev.dat                               # V(A^3) E(Ry)
    python qe_eos_fit.py ea.dat --x alat --vfactor 0.25 -o si_eos.png
    python qe_eos_fit.py --from-out scan/*/pw.out --eos murnaghan
    python qe_eos_fit.py ea.dat --x alat --xunit bohr --vfactor 0.25 --eunit Ry

Dependencies: python >= 3.9, numpy, scipy, matplotlib (Agg backend unless --show).
"""
import argparse
import re
import sys

RY_EV = 13.605693122994
BOHR_A = 0.529177210903
EVA3_GPA = 160.21766208


def bm3(v, e0, v0, b0, bp):
    eta = (v0 / v) ** (2.0 / 3.0)
    return e0 + 9.0 * v0 * b0 / 16.0 * ((eta - 1.0) ** 3 * bp + (eta - 1.0) ** 2 * (6.0 - 4.0 * eta))


def murnaghan(v, e0, v0, b0, bp):
    return e0 + b0 * v / bp * ((v0 / v) ** bp / (bp - 1.0) + 1.0) - b0 * v0 / (bp - 1.0)


def from_outputs(files):
    vs, es = [], []
    for f in files:
        v = e = None
        with open(f, errors="replace") as fh:
            for ln in fh:
                s = ln.strip()
                if s.startswith("unit-cell volume"):
                    v = float(s.split("=")[1].split()[0])
                elif s.startswith("new unit-cell volume"):
                    v = float(s.split("=")[1].split()[0])
                elif ln.startswith("!"):
                    m = re.search(r"=\s*([-+.\dEe]+)\s*Ry", ln)
                    if m:
                        e = float(m.group(1))
        if v is None or e is None:
            print("skip %s (no volume or no '!' energy)" % f, file=sys.stderr)
            continue
        vs.append(v * BOHR_A ** 3)
        es.append(e * RY_EV)
    return vs, es


def main():
    ap = argparse.ArgumentParser(description="Fit Birch-Murnaghan / Murnaghan EOS.")
    ap.add_argument("data", nargs="?", help="two-column file: V (or a) and E")
    ap.add_argument("--from-out", nargs="+", metavar="OUT", help="pw.x outputs instead of a data file")
    ap.add_argument("--x", choices=["volume", "alat"], default="volume", help="first column")
    ap.add_argument("--xunit", choices=["angstrom", "bohr"], default="angstrom")
    ap.add_argument("--vfactor", type=float, default=1.0, help="V = vfactor * a^3 (with --x alat)")
    ap.add_argument("--eunit", choices=["Ry", "eV", "Ha"], default="Ry")
    ap.add_argument("--natoms", type=int, default=1, help="divide E and V by this for per-atom output")
    ap.add_argument("--eos", choices=["bm3", "murnaghan"], default="bm3")
    ap.add_argument("-o", "--output", default="eos.png", help="figure, 'none' = no figure")
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    import numpy as np
    from scipy.optimize import curve_fit
    if args.from_out:
        V, E = from_outputs(args.from_out)
        V, E = np.array(V), np.array(E)
        a_mode = False
    else:
        if not args.data:
            ap.error("give a data file or --from-out")
        d = np.loadtxt(args.data, comments="#")
        x, E = d[:, 0], d[:, 1]
        if args.xunit == "bohr":
            x = x * BOHR_A
        a_mode = args.x == "alat"
        V = args.vfactor * x ** 3 if a_mode else x
        E = E * {"Ry": RY_EV, "eV": 1.0, "Ha": 2 * RY_EV}[args.eunit]
    if len(V) < 5:
        sys.exit("need at least 5 points (have %d)" % len(V))
    o = np.argsort(V)
    V, E = V[o] / args.natoms, E[o] / args.natoms

    c = np.polyfit(V, E, 2)
    v0 = -c[1] / (2 * c[0]) if c[0] > 0 else V[np.argmin(E)]
    b0 = max(2 * c[0] * v0, 1e-3)
    p0 = [E.min(), v0, b0, 4.0]
    func = bm3 if args.eos == "bm3" else murnaghan
    try:
        popt, pcov = curve_fit(func, V, E, p0=p0, maxfev=20000)
    except RuntimeError as exc:
        sys.exit("fit failed: %s" % exc)
    e0, v0, b0, bp = popt
    perr = np.sqrt(np.diag(pcov))
    res = E - func(V, *popt)
    print("EOS %s, %d points%s" % (args.eos, len(V), " (per atom)" if args.natoms > 1 else ""))
    print("  E0  = %.6f eV  (%.8f Ry)" % (e0, e0 / RY_EV))
    print("  V0  = %.4f A^3 (+- %.4f)  = %.4f bohr^3" % (v0, perr[1], v0 / BOHR_A ** 3))
    if a_mode:
        a0 = (v0 * args.natoms / args.vfactor) ** (1.0 / 3.0)
        print("  a0  = %.5f A = %.5f bohr" % (a0, a0 / BOHR_A))
    print("  B0  = %.2f GPa (+- %.2f)" % (b0 * EVA3_GPA, perr[2] * EVA3_GPA))
    print("  B0' = %.3f (+- %.3f)" % (bp, perr[3]))
    print("  rms residual = %.3e eV" % np.sqrt(np.mean(res ** 2)))
    if not (V.min() < v0 < V.max()):
        print("WARNING: V0 outside the data range - extend the scan", file=sys.stderr)
    elif min(v0 - V.min(), V.max() - v0) < 0.03 * v0:
        print("WARNING: minimum close to the edge of the data - add points on that side",
              file=sys.stderr)

    if args.output.lower() == "none" and not args.show:
        return
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    vv = np.linspace(V.min() * 0.98, V.max() * 1.02, 300)
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    ax.plot(V, E - e0, "o", color="#1f4e79", label="data")
    ax.plot(vv, func(vv, *popt) - e0, "-", color="#c0504d",
            label="%s: B$_0$=%.1f GPa" % (args.eos.upper(), b0 * EVA3_GPA))
    ax.axvline(v0, color="0.6", lw=0.6, ls=":")
    ax.set_xlabel("V (A$^3$%s)" % ("/atom" if args.natoms > 1 else ""))
    ax.set_ylabel("E - E$_0$ (eV%s)" % ("/atom" if args.natoms > 1 else ""))
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    if args.output.lower() != "none":
        fig.savefig(args.output, dpi=args.dpi)
        print("figure written to %s" % args.output)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
