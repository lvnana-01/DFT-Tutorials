#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qe_workfunction.py -- work function of a slab from average.x output (avg.dat) or a cube file.
                      读取 average.x 的 avg.dat 或 pp.x 写出的 cube（plot_num=11, output_format=6），
                      求平面平均势、真空能级 V_vac 与功函数 Phi = V_vac - E_F，并作图。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Workflow (chapter 9)
  1. pw.x scf of the slab (vacuum >= 15 A, dipole correction for asymmetric slabs:
     tefield=.true., dipfield=.true., edir, emaxpos, eopreg, eamp=0).
  2. pp.x  &INPUTPP plot_num=11 (V_bare + V_H) filplot='slab.pot' /
           &PLOT iflag=3, output_format=6, fileout='slab_pot.cube' /   (cube, optional)
  3. average.x  < avg.in   (reads filplot, writes avg.dat)
  4. python qe_workfunction.py avg.dat --scf-out slab.scf.out
avg.dat format (PP/src/average.f90): three columns "z  p(z)  m(z)": z in bohr (from 0 to the
  cell length along idir), planar average p and macroscopic average m in the units of the
  pp.x quantity -> Ry for plot_num=11.  QE potentials are electron potential energies, so the
  vacuum plateau lies above E_F and Phi = V_vac - E_F > 0.
Cube (Gaussian format written by pp.x): lengths in bohr, values in Ry for plot_num=11; the
  planar average is taken over the plane perpendicular to --axis.
Vacuum level: average of p(z) inside --window Z1 Z2 (angstrom), or automatically in the
  flattest region: the central 20 % of the widest region where |dp/dz| < --flat-tol eV/A
  (two plateaux are reported when a dipole correction creates a step).

Usage examples
--------------
    python qe_workfunction.py avg.dat --scf-out pt111.scf.out -o pt111_wf.png
    python qe_workfunction.py avg.dat --efermi 5.123 --window 22 26
    python qe_workfunction.py slab_pot.cube --axis 3 --scf-out slab.scf.out --csv vz.csv
    python qe_workfunction.py avg.dat --efermi 0.0 --column macro     # use m(z)

Dependencies: python >= 3.9, numpy, matplotlib (Agg backend unless --show).
"""
import argparse
import re
import sys

RY_EV = 13.605693122994
BOHR_A = 0.529177210903
NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"


def read_cube_planar(path, axis):
    import numpy as np
    with open(path) as fh:
        fh.readline()
        fh.readline()
        v = fh.readline().split()
        nat = int(v[0])
        n, vec = [], []
        for _ in range(3):
            v = fh.readline().split()
            n.append(int(v[0]))
            vec.append([float(x) for x in v[1:4]])
        if any(k < 0 for k in n):
            sys.exit("cube lengths in angstrom (negative N) are not expected from pp.x")
        for _ in range(abs(nat)):
            fh.readline()
        data = np.array(fh.read().split(), dtype=float)
    if data.size != n[0] * n[1] * n[2]:
        sys.exit("cube data size %d != %d x %d x %d" % (data.size, n[0], n[1], n[2]))
    data = data.reshape(n)
    ax = axis - 1
    others = tuple(i for i in range(3) if i != ax)
    prof = data.mean(axis=others)
    length = np.linalg.norm(vec[ax]) * n[ax]  # bohr
    z = np.arange(n[ax]) * length / n[ax]
    return z, prof, nat


def ef_from_scf(path):
    ef = None
    with open(path, errors="replace") as fh:
        for ln in fh:
            s = ln.strip()
            if s.startswith("the Fermi energy is"):
                ef = float(s.split()[4])
            elif s.startswith("highest occupied, lowest unoccupied level") or \
                    s.startswith("highest occupied level"):
                ef = float(re.findall(NUM, s.split(":")[1])[0])
    if ef is None:
        sys.exit("no Fermi energy / HOMO in %s" % path)
    return ef


def plateaux(z, v, tol):
    """Return list of (z_start, z_end, mean) flat regions, widest first."""
    import numpy as np
    dv = np.abs(np.gradient(v, z))
    flat = dv < tol
    regs, start = [], None
    for i, f in enumerate(np.append(flat, False)):
        if f and start is None:
            start = i
        elif not f and start is not None:
            if i - start >= 3:
                regs.append((start, i - 1))
            start = None
    # periodic wrap: join first and last region
    if len(regs) > 1 and regs[0][0] == 0 and regs[-1][1] == len(z) - 1:
        a, b = regs.pop(-1), regs.pop(0)
        regs.append((a[0], b[1] + len(z)))
    out = []
    for a, b in regs:
        idx = np.arange(a, b + 1) % len(z)
        w = b - a + 1
        core = idx[int(0.4 * w):max(int(0.6 * w), int(0.4 * w) + 1)]
        out.append((z[idx[0]], z[idx[-1]], float(np.mean(v[core])), w))
    out.sort(key=lambda r: -r[3])
    return out


def main():
    ap = argparse.ArgumentParser(description="Work function from average.x avg.dat or a cube file.")
    ap.add_argument("data", help="avg.dat (average.x) or .cube (pp.x output_format=6)")
    ap.add_argument("--axis", type=int, choices=[1, 2, 3], default=3, help="cube: surface normal")
    ap.add_argument("--column", choices=["planar", "macro"], default="planar",
                    help="avg.dat: use p(z) (col 2) or macroscopic average m(z) (col 3)")
    ap.add_argument("--efermi", type=float, help="E_F in eV")
    ap.add_argument("--scf-out", help="pw.x scf output (E_F)")
    ap.add_argument("--window", type=float, nargs=2, metavar=("Z1", "Z2"),
                    help="vacuum window in angstrom (else automatic)")
    ap.add_argument("--flat-tol", type=float, default=0.05, help="eV/A flatness for auto vacuum")
    ap.add_argument("--csv", help="write z(A), V(eV) to CSV")
    ap.add_argument("-o", "--output", default="workfunction.png", help="figure, 'none' = no figure")
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    import numpy as np
    if args.data.lower().endswith(".cube"):
        z, v, _ = read_cube_planar(args.data, args.axis)
    else:
        d = np.loadtxt(args.data)
        z, v = d[:, 0], d[:, 1 if args.column == "planar" else 2]
    z = z * BOHR_A
    v = v * RY_EV
    if args.efermi is not None:
        ef = args.efermi
    elif args.scf_out:
        ef = ef_from_scf(args.scf_out)
    else:
        ef = None
        print("WARNING: no E_F (--efermi/--scf-out): only the vacuum level is reported",
              file=sys.stderr)
    print("profile: %d points, L = %.3f A, V from %.3f to %.3f eV" % (len(z), z[-1] + z[1] - z[0],
                                                                       v.min(), v.max()))
    levels = []
    if args.window:
        z1, z2 = sorted(args.window)
        m = (z >= z1) & (z <= z2)
        if not m.any():
            sys.exit("window contains no points")
        vv = float(np.mean(v[m]))
        print("vacuum window %.2f-%.2f A: V_vac = %.4f eV (spread %.4f eV)"
              % (z1, z2, vv, v[m].max() - v[m].min()))
        levels.append((z1, z2, vv))
    else:
        regs = plateaux(z, v, args.flat_tol)
        if not regs:
            sys.exit("no flat vacuum region found; increase --flat-tol or give --window")
        top = [r for r in regs if r[2] > v.mean()][:2] or regs[:1]
        for a, b, vv, _ in top:
            print("vacuum plateau %.2f-%.2f A: V_vac = %.4f eV" % (a, b, vv))
            levels.append((a, b, vv))
        if len(top) == 2 and abs(top[0][2] - top[1][2]) > 0.05:
            print("two vacuum levels differ by %.3f eV (dipole correction / asymmetric slab)"
                  % abs(top[0][2] - top[1][2]))
    if ef is not None:
        print("E_F = %.4f eV" % ef)
        for a, b, vv in levels:
            print("work function Phi = V_vac - E_F = %.4f eV  (vacuum %.1f-%.1f A)" % (vv - ef, a, b))
    if args.csv:
        np.savetxt(args.csv, np.column_stack([z, v]), delimiter=",", fmt="%.6f",
                   header="z_A,V_eV", comments="")
        print("CSV written to %s" % args.csv)
    if args.output.lower() == "none" and not args.show:
        return
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    ax.plot(z, v, color="#1f4e79", lw=1.1, label="planar average" if args.column == "planar"
            else "macroscopic average")
    for a, b, vv in levels:
        ax.hlines(vv, a, b, color="#c0504d", lw=2)
        ax.text((a + b) / 2, vv, " V$_{vac}$=%.2f" % vv, va="bottom", ha="center", fontsize=8)
    if ef is not None:
        ax.axhline(ef, color="0.3", ls="--", lw=0.8, label="E$_F$")
    ax.set_xlabel("z (A)")
    ax.set_ylabel("V (eV)")
    ax.set_xlim(z[0], z[-1])
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    if args.output.lower() != "none":
        fig.savefig(args.output, dpi=args.dpi)
        print("figure written to %s" % args.output)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
