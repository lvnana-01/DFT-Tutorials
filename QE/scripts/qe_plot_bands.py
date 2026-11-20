#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qe_plot_bands.py -- plot QE band structures from bands.x (filband.gnu) or the XML file.
                    读取 bands.x 输出的 <filband>.gnu 或 pw.x 的 data-file-schema.xml，
                    以费米能/VBM 为零点，按 crystal_b 标签画能带，并给出带隙（直接/间接）。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Input formats
-------------
  <filband>.gnu   written by bands.x: two columns "k  E(eV)", one block per band separated
                  by a blank line; k = cumulative path length in units of 2*pi/alat.
  bands.x stdout  lines "high-symmetry point:  kx ky kz   x coordinate   X"  -> tick positions
                  (give with --bandsx-out; a jump in the path appears as two points with the same x).
  data-file-schema.xml  (<outdir>/<prefix>.save/) of the 'bands' run: <ks_energies> blocks with
                  <k_point> (cartesian, 2*pi/alat) and <eigenvalues> (hartree); <lsda> true ->
                  first nbnd_up values are spin up, then nbnd_dw spin down.
Energy zero (in this order): --efermi EV, value from --scf-out (pw.x scf/nscf output:
  "the Fermi energy is" or "highest occupied, lowest unoccupied level (ev):" /
  "highest occupied level (ev):"), or the XML <fermi_energy>/<highestOccupiedLevel>.
  Do NOT take E_F from the 'bands' run itself (it is not self-consistent).
Gap: VBM = max(E <= E_ref + tol), CBM = min(E > E_ref + tol) over all plotted k; reported
  as direct/indirect with the k index/x position.  For metals the "gap" is meaningless.

Usage examples
--------------
    python qe_plot_bands.py si.bands.dat.gnu --bandsx-out si.bands_pp.out --scf-out si.scf.out \
           --labels G X W K G L --ylim -12 6 -o si_bands.png
    python qe_plot_bands.py si.bands.dat.gnu --efermi 6.2 --labels "G X|U K G L" --csv bands.csv
    python qe_plot_bands.py tmp/si.save/data-file-schema.xml --labels G X W K G L
    python qe_plot_bands.py fe_up.gnu --spin2 fe_dw.gnu --scf-out fe.scf.out   # spin_component 1/2

Labels: one per tick (high-symmetry point).  "X|U" makes a single tick labelled X|U at a path
jump; "G" is printed as Gamma.
Dependencies: python >= 3.9, numpy, matplotlib (Agg backend unless --show).
"""
import argparse
import re
import sys

HA_EV = 27.211386245988
NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"


def read_gnu(path):
    import numpy as np
    bands, cur = [], []
    with open(path) as fh:
        for ln in fh:
            v = ln.split()
            if len(v) < 2:
                if cur:
                    bands.append(cur)
                    cur = []
                continue
            cur.append((float(v[0]), float(v[1])))
    if cur:
        bands.append(cur)
    if not bands:
        sys.exit("no data in %s" % path)
    nk = min(len(b) for b in bands)
    arr = np.array([[p for p in b[:nk]] for b in bands])  # (nbnd, nk, 2)
    return arr[0, :, 0], arr[:, :, 1].T  # x(nk), E(nk, nbnd)


def read_ticks(path):
    xs = []
    with open(path, errors="replace") as fh:
        for ln in fh:
            if "high-symmetry point:" in ln and "x coordinate" in ln:
                xs.append(float(ln.split("x coordinate")[1].split()[0]))
    return xs


def read_xml(path):
    import numpy as np
    import xml.etree.ElementTree as ET
    root = ET.parse(path).getroot()
    bs = root.find(".//band_structure")
    if bs is None:
        sys.exit("no <band_structure> in %s" % path)
    lsda = (bs.findtext("lsda") or "false").strip().lower() == "true"
    ks = bs.findall("ks_energies")
    kpts = np.array([[float(x) for x in k.find("k_point").text.split()] for k in ks])
    eig = np.array([[float(x) for x in k.find("eigenvalues").text.split()] for k in ks]) * HA_EV
    ef = None
    for tag in ("fermi_energy", "highestOccupiedLevel"):
        t = bs.findtext(tag)
        if t:
            ef = float(t) * HA_EV
            break
    dk = np.linalg.norm(np.diff(kpts, axis=0), axis=1)
    x = np.concatenate([[0.0], np.cumsum(dk)])
    if lsda:
        nup = int(bs.findtext("nbnd_up"))
        return x, eig[:, :nup], eig[:, nup:], ef, kpts
    return x, eig, None, ef, kpts


def ticks_from_kpts(x, kpts):
    """High-symmetry points from direction changes of the k path (XML input)."""
    import numpy as np
    t = [x[0]]
    d = np.diff(kpts, axis=0)
    for i in range(1, len(d)):
        a, b = d[i - 1], d[i]
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        if na < 1e-8 or nb < 1e-8:
            continue
        if abs(np.dot(a, b) / (na * nb) - 1.0) > 1e-4 or abs(nb / na - 1) > 0.5:
            t.append(x[i])
    t.append(x[-1])
    return t


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
            elif s.startswith("highest occupied, lowest unoccupied level"):
                ef = float(re.findall(NUM, s.split(":")[1])[0])
            elif s.startswith("highest occupied level"):
                ef = float(re.findall(NUM, s.split(":")[1])[0])
    if ef is None:
        sys.exit("no Fermi energy / HOMO found in %s" % path)
    return ef


def gap(x, e_list, tol):
    import numpy as np
    allE = np.concatenate(e_list, axis=1)
    occ = np.where(allE <= tol, allE, -np.inf)
    emp = np.where(allE > tol, allE, np.inf)
    vk = occ.max(axis=1)
    ck = emp.min(axis=1)
    iv, ic = int(np.argmax(vk)), int(np.argmin(ck))
    if not np.isfinite(vk[iv]) or not np.isfinite(ck[ic]):
        return None
    g = ck[ic] - vk[iv]
    gd = ck - vk
    idir = int(np.argmin(gd))
    # metal check: a band crossing zero
    crossing = any((b.min() < -tol) and (b.max() > tol) for b in allE.T)
    return dict(vbm=vk[iv], cbm=ck[ic], gap=g, ik_v=iv, ik_c=ic, x_v=x[iv], x_c=x[ic],
                direct_gap=gd[idir], ik_direct=idir, metal=crossing)


def main():
    ap = argparse.ArgumentParser(description="Plot QE band structure (bands.x .gnu or XML).")
    ap.add_argument("data", help="<filband>.gnu from bands.x or data-file-schema.xml")
    ap.add_argument("--spin2", help="second .gnu file (spin down, bands.x spin_component=2)")
    ap.add_argument("--bandsx-out", help="bands.x stdout (for high-symmetry x coordinates)")
    ap.add_argument("--ticks", type=float, nargs="+", help="explicit tick x positions (2pi/alat)")
    ap.add_argument("--labels", nargs="+",
                    help="tick labels, e.g. G X W K G L (one string with spaces also accepted)")
    ap.add_argument("--efermi", type=float, help="reference energy in eV")
    ap.add_argument("--scf-out", help="pw.x scf/nscf output to read E_F or HOMO from")
    ap.add_argument("--no-shift", action="store_true", help="do not shift energies")
    ap.add_argument("--tol", type=float, default=1e-3, help="eV tolerance for occupied/empty")
    ap.add_argument("--ylim", type=float, nargs=2, default=[-6.0, 6.0])
    ap.add_argument("--title")
    ap.add_argument("--csv", help="write x and shifted energies to CSV")
    ap.add_argument("-o", "--output", default="bands.png", help="figure file, 'none' = no figure")
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--show", action="store_true", help="interactive window")
    args = ap.parse_args()

    import numpy as np
    e_dn = None
    ef_file = None
    if args.data.lower().endswith(".xml"):
        x, e_up, e_dn, ef_file, kpts = read_xml(args.data)
        ticks = args.ticks or ticks_from_kpts(x, kpts)
    else:
        x, e_up = read_gnu(args.data)
        if args.spin2:
            x2, e_dn = read_gnu(args.spin2)
            if len(x2) != len(x):
                sys.exit("spin up/down files have different k points")
        if args.ticks:
            ticks = args.ticks
        elif args.bandsx_out:
            ticks = read_ticks(args.bandsx_out)
        else:
            ticks = [x[0], x[-1]]
            print("NOTE: no --bandsx-out/--ticks: only the path ends are marked", file=sys.stderr)
    # merge duplicate tick positions (path jumps)
    merged = []
    for t in ticks:
        if not merged or abs(t - merged[-1]) > 1e-4:
            merged.append(t)
    ticks = merged

    if args.efermi is not None:
        ef = args.efermi
    elif args.scf_out:
        ef = ef_from_scf(args.scf_out)
    elif ef_file is not None:
        ef = ef_file
        print("NOTE: using E_F/HOMO from the XML file (%.4f eV); for a 'bands' run prefer "
              "--scf-out of the scf run" % ef, file=sys.stderr)
    else:
        ef = 0.0
        print("WARNING: no reference energy given; energies are absolute", file=sys.stderr)
    shift = 0.0 if args.no_shift else ef
    sets = [e_up - shift] + ([e_dn - shift] if e_dn is not None else [])

    print("k points %d, bands %d%s, E_ref = %.4f eV" % (len(x), e_up.shape[1],
                                                       " (x2 spin)" if e_dn is not None else "", ef))
    g = gap(x, sets, args.tol + (0.0 if not args.no_shift else ef))
    if g:
        if g["metal"]:
            print("bands cross E_ref -> metallic (no gap)")
        else:
            kind = "direct" if g["ik_v"] == g["ik_c"] else "indirect"
            print("VBM %.4f eV at k#%d (x=%.4f), CBM %.4f eV at k#%d (x=%.4f)"
                  % (g["vbm"], g["ik_v"] + 1, g["x_v"], g["cbm"], g["ik_c"] + 1, g["x_c"]))
            print("band gap %.4f eV (%s); smallest direct gap %.4f eV at k#%d"
                  % (g["gap"], kind, g["direct_gap"], g["ik_direct"] + 1))
            print("NOTE: gap on the plotted path only; DFT (LDA/GGA) underestimates gaps.")

    if args.csv:
        with open(args.csv, "w") as fh:
            cols = ["x"] + ["up%d" % (i + 1) if e_dn is not None else "b%d" % (i + 1)
                            for i in range(e_up.shape[1])]
            if e_dn is not None:
                cols += ["dn%d" % (i + 1) for i in range(e_dn.shape[1])]
            fh.write(",".join(cols) + "\n")
            data = np.column_stack([x] + sets)
            for row in data:
                fh.write(",".join("%.6f" % v for v in row) + "\n")
        print("CSV written to %s" % args.csv)

    if args.output.lower() == "none" and not args.show:
        return
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    ax.plot(x, sets[0], color="#1f4e79", lw=1.2)
    if e_dn is not None:
        ax.plot(x, sets[1], color="#c0504d", lw=1.0, ls="--")
        ax.plot([], [], color="#1f4e79", label="spin up")
        ax.plot([], [], color="#c0504d", ls="--", label="spin down")
        ax.legend(frameon=False, fontsize=8)
    for t in ticks:
        ax.axvline(t, color="0.6", lw=0.6)
    if not args.no_shift:
        ax.axhline(0.0, color="0.3", lw=0.6, ls=":")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(*args.ylim)
    if args.labels:
        labels = args.labels[0].split() if len(args.labels) == 1 else args.labels
        labels = [r"$\Gamma$" if lab.upper() in ("G", "GAMMA") else lab for lab in labels]
        if len(labels) == len(ticks):
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
        else:
            print("WARNING: %d labels but %d ticks %s - labels ignored"
                  % (len(labels), len(ticks), ["%.4f" % t for t in ticks]), file=sys.stderr)
            ax.set_xticks(ticks)
            ax.set_xticklabels([])
    else:
        ax.set_xticks(ticks)
        ax.set_xticklabels([])
    ax.set_ylabel("E - E$_{ref}$ (eV)" if not args.no_shift else "E (eV)")
    if args.title:
        ax.set_title(args.title)
    fig.tight_layout()
    if args.output.lower() != "none":
        fig.savefig(args.output, dpi=args.dpi)
        print("figure written to %s" % args.output)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
