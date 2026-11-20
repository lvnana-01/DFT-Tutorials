#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qe_plot_dos.py -- plot total DOS (dos.x) and projected DOS (projwfc.x) of Quantum ESPRESSO.
                  读取 dos.x 的 <fildos> 与 projwfc.x 的 <filpdos>.pdos_tot / .pdos_atm#N(El)_wfc#M(l)
                  文件，按元素 / 元素+轨道 / 原子求和，自旋分辨作图（自旋向下画成负值）。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

File formats (QE 7.x, PP/src/dos.f90 and PP/src/partialdos.f90)
  dos.x       "#  E (eV)   dos(E)     Int dos(E) EFermi =    6.123 eV"
              columns E, dos, int  (nspin=1)  or  E, dosup, dosdw, int  (nspin=2)
  pdos_tot    "# E (eV)  dos(E)    pdos(E)"            (nspin=1, also noncollinear)
              "# E (eV)  dosup(E)   dosdw(E)  pdosup(E)  pdosdw(E)"   (nspin=2)
  pdos_atm#1(Si)_wfc#2(p)
              "# E (eV)   ldos(E)  pdos(E) pdos(E) pdos(E)"           (nspin=1)
              "# E (eV)  ldosup(E)  ldosdw(E) pdosup(E) pdosdw(E) ..."   (nspin=2)
              with spin-orbit the name ends in "(p_j1.5)" and only ldos has meaning per j.
  The first data column after E is the m-summed (l)dos, which is what this script sums.
  k-resolved files (kresolveddos=.true., extra "ik" column) are not supported.
E_F: from --efermi, the "EFermi =" field of the dos.x header, or --scf-out
     (pw.x nscf output: "the Fermi energy is" / "highest occupied, lowest unoccupied level").

Usage examples
--------------
    python qe_plot_dos.py --dos si.dos --xlim -12 6 -o si_dos.png
    python qe_plot_dos.py --pdos-prefix si.pdos --by element --scf-out si.nscf.out
    python qe_plot_dos.py --pdos-prefix fe.pdos --by orbital --efermi 17.8 --xlim -8 4
    python qe_plot_dos.py --dos tio2.dos --pdos-prefix tio2.pdos --by orbital \
           --select Ti-d O-p --csv pdos.csv
    python qe_plot_dos.py --pdos-prefix nio.pdos --by atom --select 1 2   # atoms #1, #2

--by element : sum over all atoms and l of each element   (labels "Fe")
--by orbital : sum over atoms of each element, per l       (labels "Fe-d")
--by atom    : per atom, all l                             (labels "1-Fe")
Dependencies: python >= 3.9, numpy, matplotlib (Agg backend unless --show).
"""
import argparse
import glob
import os
import re
import sys

NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
PAT = re.compile(r"\.pdos_atm#(\d+)\(([^)]+)\)_wfc#(\d+)\(([spdf])(?:_j([\d.]+))?\)$")


def load(path):
    import numpy as np
    header = ""
    with open(path) as fh:
        for ln in fh:
            if ln.lstrip().startswith("#"):
                header = ln
                break
    if "ik" in header.split():
        sys.exit("%s is k-resolved (kresolveddos) - not supported" % path)
    data = np.loadtxt(path, comments="#")
    return header, data


def read_dos(path):
    header, d = load(path)
    ef = None
    m = re.search(r"EFermi\s*=\s*(" + NUM + r")(?:\s+(" + NUM + r"))?", header)
    if m:
        ef = float(m.group(1)) if m.group(2) is None else max(float(m.group(1)), float(m.group(2)))
    spin = "dosup" in header
    if spin:
        return d[:, 0], d[:, 1], d[:, 2], ef
    return d[:, 0], d[:, 1], None, ef


def read_pdos(prefix):
    files = sorted(glob.glob(glob.escape(prefix) + ".pdos_atm#*"))
    if not files:
        sys.exit("no files %s.pdos_atm#* found" % prefix)
    out = []
    for f in files:
        m = PAT.search(f)
        if not m:
            print("skip %s (name not recognised)" % f, file=sys.stderr)
            continue
        header, d = load(f)
        spin = "ldosup" in header
        rec = dict(atom=int(m.group(1)), el=m.group(2), wfc=int(m.group(3)), l=m.group(4),
                   j=m.group(5), E=d[:, 0], up=d[:, 1], dn=d[:, 2] if spin else None)
        out.append(rec)
    return out


def ef_from_scf(path):
    ef = None
    with open(path, errors="replace") as fh:
        for ln in fh:
            s = ln.strip()
            if s.startswith("the Fermi energy is"):
                ef = float(s.split()[4])
            elif s.startswith("the spin up/dw Fermi energies are"):
                v = re.findall(NUM, s)
                ef = max(float(v[0]), float(v[1]))
            elif s.startswith("highest occupied, lowest unoccupied level") or \
                    s.startswith("highest occupied level"):
                ef = float(re.findall(NUM, s.split(":")[1])[0])
    if ef is None:
        sys.exit("no Fermi energy / HOMO in %s" % path)
    return ef


def group(recs, by):
    import numpy as np
    groups = {}
    for r in recs:
        if by == "element":
            key = r["el"]
        elif by == "orbital":
            key = "%s-%s" % (r["el"], r["l"])
        else:
            key = "%d-%s" % (r["atom"], r["el"])
        g = groups.setdefault(key, {"E": r["E"], "up": np.zeros_like(r["up"]),
                                    "dn": None if r["dn"] is None else np.zeros_like(r["dn"])})
        if len(r["E"]) != len(g["E"]):
            sys.exit("PDOS files have different energy grids")
        g["up"] += r["up"]
        if r["dn"] is not None:
            g["dn"] += r["dn"]
    return groups


def main():
    ap = argparse.ArgumentParser(description="Plot DOS (dos.x) and PDOS (projwfc.x).")
    ap.add_argument("--dos", help="dos.x output file (fildos)")
    ap.add_argument("--pdos-prefix", help="filpdos of projwfc.x (files <prefix>.pdos_atm#...)")
    ap.add_argument("--by", choices=["element", "orbital", "atom"], default="element")
    ap.add_argument("--select", nargs="+", help="only these groups (e.g. Fe-d O-p, or atom numbers)")
    ap.add_argument("--pdos-tot", action="store_true",
                    help="also plot the total from <prefix>.pdos_tot (sum of projections)")
    ap.add_argument("--efermi", type=float, help="E_F in eV (overrides file values)")
    ap.add_argument("--scf-out", help="pw.x (n)scf output to read E_F from")
    ap.add_argument("--no-shift", action="store_true", help="absolute energies")
    ap.add_argument("--xlim", type=float, nargs=2, default=[-10.0, 5.0])
    ap.add_argument("--ymax", type=float, help="upper DOS limit")
    ap.add_argument("--csv", help="write all plotted curves to CSV")
    ap.add_argument("-o", "--output", default="dos.png", help="figure file, 'none' = no figure")
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()
    if not args.dos and not args.pdos_prefix:
        ap.error("give --dos and/or --pdos-prefix")

    import numpy as np
    curves = []  # (label, E, up, dn, style)
    ef_file = None
    if args.dos:
        E, up, dn, ef_file = read_dos(args.dos)
        curves.append(("total (dos.x)", E, up, dn, dict(color="k", lw=1.2)))
    if args.pdos_prefix:
        recs = read_pdos(args.pdos_prefix)
        if any(r["j"] for r in recs):
            print("NOTE: spin-orbit PDOS: projections of the same l and different j are summed",
                  file=sys.stderr)
        groups = group(recs, args.by)
        keys = list(groups)
        if args.select:
            if args.by == "atom":
                keys = [k for k in keys if k.split("-")[0] in args.select]
            else:
                keys = [k for k in keys if k in args.select]
            if not keys:
                sys.exit("--select matched nothing; available: %s" % ", ".join(groups))
        for k in keys:
            g = groups[k]
            curves.append((k, g["E"], g["up"], g["dn"], dict(lw=1.0)))
        tot = args.pdos_prefix + ".pdos_tot"
        if args.pdos_tot and os.path.isfile(tot):
            h, d = load(tot)
            if "dosup" in h.replace("pdosup", ""):   # nspin=2: E dosup dosdw pdosup pdosdw
                curves.append(("sum pdos", d[:, 0], d[:, 3], d[:, 4], dict(color="0.5", ls="--")))
            elif "pdosup" in h:                       # noncollinear: E dos pdosup pdosdw
                curves.append(("sum pdos", d[:, 0], d[:, 2], d[:, 3], dict(color="0.5", ls="--")))
            else:
                curves.append(("sum pdos", d[:, 0], d[:, 2], None, dict(color="0.5", ls="--")))
        print("PDOS groups (%s): %s" % (args.by, ", ".join(groups)))

    if args.efermi is not None:
        ef = args.efermi
    elif args.scf_out:
        ef = ef_from_scf(args.scf_out)
    elif ef_file is not None:
        ef = ef_file
    else:
        ef = 0.0
        if not args.no_shift:
            print("WARNING: no E_F (use --efermi/--scf-out/--dos); energies not shifted",
                  file=sys.stderr)
    shift = 0.0 if args.no_shift else ef
    print("E_F = %.4f eV (%s)" % (ef, "shifted to 0" if shift else "not shifted"))

    spin = any(c[3] is not None for c in curves)
    for lab, E, up, dn, _ in curves:
        win = (E - shift >= args.xlim[0]) & (E - shift <= args.xlim[1])
        de = np.gradient(E)
        if spin and dn is not None:
            print("  %-16s integral in window: up %.3f  down %.3f" % (lab, np.sum(up[win] * de[win]),
                                                                    np.sum(dn[win] * de[win])))
        else:
            print("  %-16s integral in window: %.3f" % (lab, np.sum(up[win] * de[win])))

    if args.csv:
        E0 = curves[0][1]
        cols, data = ["E-EF"], [E0 - shift]
        for lab, E, up, dn, _ in curves:
            if len(E) != len(E0):
                print("CSV: %s skipped (other energy grid)" % lab, file=sys.stderr)
                continue
            cols.append(lab + ("_up" if dn is not None else ""))
            data.append(up)
            if dn is not None:
                cols.append(lab + "_dn")
                data.append(dn)
        np.savetxt(args.csv, np.column_stack(data), delimiter=",", header=",".join(cols),
                   comments="", fmt="%.6e")
        print("CSV written to %s" % args.csv)

    if args.output.lower() == "none" and not args.show:
        return
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for lab, E, up, dn, st in curves:
        line, = ax.plot(E - shift, up, label=lab, **st)
        if dn is not None:
            ax.plot(E - shift, -dn, color=line.get_color(), **{k: v for k, v in st.items()
                                                               if k != "color"})
    if spin:
        ax.axhline(0, color="0.3", lw=0.5)
    if not args.no_shift:
        ax.axvline(0, color="0.3", lw=0.6, ls=":")
    ax.set_xlim(*args.xlim)
    if args.ymax:
        ax.set_ylim(-args.ymax if spin else 0, args.ymax)
    elif not spin:
        ax.set_ylim(bottom=0)
    ax.set_xlabel("E - E$_F$ (eV)" if not args.no_shift else "E (eV)")
    ax.set_ylabel("DOS (states/eV)" + ("  (down < 0)" if spin else ""))
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    if args.output.lower() != "none":
        fig.savefig(args.output, dpi=args.dpi)
        print("figure written to %s" % args.output)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
