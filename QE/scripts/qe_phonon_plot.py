#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qe_phonon_plot.py -- phonon dispersion and DOS from matdyn.x output; flags imaginary modes.
                     读取 matdyn.x 的 <flfrq>.gp（色散）与 <fldos>（声子 DOS），画色散 + DOS，
                     标出虚频（matdyn.x 把 omega^2 < 0 的模式写成负频率）。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

File formats (PHonon/PH/matdyn.f90, QE 7.x)
  <flfrq>.gp   "pathL  w_1 ... w_3nat" per q point; frequencies in cm^-1, negative = imaginary;
               pathL = cumulative |q_n - q_(n-1)| in 2*pi/alat (no break at path jumps).
  <flfrq>      " &plot nbnd= 6, nks= 151 /" then q (cartesian, 2*pi/alat) and frequencies.
  <fldos>      "# Frequency[cm^-1] DOS PDOS" then E, total DOS, per-atom PDOS (dos=.true.).
  matdyn.x input with q_in_band_form=.true.: after the namelist, nq and nq lines
               "qx qy qz  n  [! label]"  -> tick positions are the cumulative point indices.
Ticks: from --matdyn-in (recommended; labels taken from "! X" comments if --labels is not
       given), or from direction changes of the q list in <flfrq> (--freq-file).

Usage examples
--------------
    python qe_phonon_plot.py --gp si.freq.gp --matdyn-in matdyn.in --dos si.dos -o si_ph.png
    python qe_phonon_plot.py --gp si.freq.gp --freq-file si.freq --labels G X U K G L
    python qe_phonon_plot.py --gp mos2.freq.gp --matdyn-in matdyn.in --unit THz --imag-tol 5
    python qe_phonon_plot.py --dos si.dos --unit meV                # DOS only

Units: cm^-1 (default), THz (1 THz = 33.35641 cm^-1), meV (1 meV = 8.065544 cm^-1).
Imaginary modes: any frequency < -imag_tol (cm^-1) is reported with its q index; small
negative acoustic values near Gamma usually mean ASR/convergence issues (asr='crystal',
tighter tr2_ph, denser q grid), large ones a dynamical instability.
Dependencies: python >= 3.9, numpy, matplotlib (Agg backend unless --show).
"""
import argparse
import re
import sys

CONV = {"cm-1": 1.0, "THz": 1.0 / 33.35641, "meV": 1.0 / 8.0655439}
NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"


def read_matdyn_in(path):
    """Return (counts, labels) of a q_in_band_form list, or (None, None)."""
    with open(path, errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    k = 0
    while k < len(lines) and lines[k].split("!")[0].strip() not in ("/",) \
            and not lines[k].split("!")[0].strip().endswith("/"):
        k += 1
    k += 1
    rows = [ln for ln in lines[k:] if ln.split("!")[0].strip()]
    if not rows:
        return None, None
    try:
        nq = int(rows[0].split()[0])
    except ValueError:
        return None, None
    counts, labels = [], []
    for ln in rows[1:1 + nq]:
        body, _, com = ln.partition("!")
        v = body.split()
        counts.append(int(float(v[3])) if len(v) > 3 else 1)
        labels.append(com.strip().split()[0] if com.strip() else "")
    return counts, labels


def read_freq_q(path):
    import numpy as np
    with open(path) as fh:
        txt = fh.read()
    m = re.search(r"nbnd=\s*(\d+)\s*,\s*nks=\s*(\d+)", txt)
    if not m:
        sys.exit("%s: no '&plot nbnd=..., nks=...' header" % path)
    nb, nk = int(m.group(1)), int(m.group(2))
    nums = [float(x) for x in re.findall(NUM, txt[txt.index("/", m.end()) + 1:])]
    per = (nb + 4) if len(nums) == nk * (nb + 4) else (nb + 3)
    q = np.array([nums[i * per:i * per + 3] for i in range(nk)])
    return q


def main():
    ap = argparse.ArgumentParser(description="Plot phonon dispersion/DOS from matdyn.x.")
    ap.add_argument("--gp", help="<flfrq>.gp file (dispersion)")
    ap.add_argument("--dos", help="<fldos> file (phonon DOS, matdyn.x dos=.true.)")
    ap.add_argument("--matdyn-in", help="matdyn.x input with q_in_band_form list (for ticks)")
    ap.add_argument("--freq-file", help="<flfrq> file (q list) to locate high-symmetry points")
    ap.add_argument("--labels", nargs="+", help="tick labels, G = Gamma")
    ap.add_argument("--unit", choices=list(CONV), default="cm-1")
    ap.add_argument("--imag-tol", type=float, default=1.0, help="cm^-1; |w| below = numerical noise")
    ap.add_argument("--ylim", type=float, nargs=2)
    ap.add_argument("--pdos", action="store_true", help="also plot per-atom PDOS columns")
    ap.add_argument("--csv", help="write the dispersion (converted units) to CSV")
    ap.add_argument("-o", "--output", default="phonon.png", help="figure, 'none' = no figure")
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()
    if not args.gp and not args.dos:
        ap.error("give --gp and/or --dos")

    import numpy as np
    f = CONV[args.unit]
    x = w = None
    ticks = []
    labels = args.labels[0].split() if args.labels and len(args.labels) == 1 else args.labels
    if args.gp:
        d = np.loadtxt(args.gp)
        if d.ndim == 1:
            d = d[None, :]
        x, w = d[:, 0], d[:, 1:]
        print("dispersion: %d q points, %d modes" % (len(x), w.shape[1]))
        if args.matdyn_in:
            counts, lab_in = read_matdyn_in(args.matdyn_in)
            if counts:
                idx = np.concatenate([[0], np.cumsum(counts[:-1])]).astype(int)
                idx = idx[idx < len(x)]
                ticks = list(x[idx])
                if not labels and any(lab_in):
                    labels = lab_in[:len(ticks)]
                # NOTE: matdyn.x does not break the .gp path at jumps (segment with n=1):
                # the jump appears as a short segment between two ticks.
        elif args.freq_file:
            q = read_freq_q(args.freq_file)
            dq = np.diff(q, axis=0)
            ticks = [x[0]]
            for i in range(1, len(dq)):
                a, b = dq[i - 1], dq[i]
                na, nb = np.linalg.norm(a), np.linalg.norm(b)
                if na > 1e-8 and nb > 1e-8 and abs(np.dot(a, b) / (na * nb) - 1) > 1e-4:
                    ticks.append(x[i])
            ticks.append(x[-1])
        else:
            ticks = [x[0], x[-1]]
        neg = np.argwhere(w < -abs(args.imag_tol))
        if len(neg):
            qs = sorted(set(int(i) for i in neg[:, 0]))
            print("IMAGINARY modes (< -%.1f cm^-1): min %.2f cm^-1, at %d q point(s): %s%s"
                  % (args.imag_tol, w.min(), len(qs), ", ".join(str(i + 1) for i in qs[:15]),
                     " ..." if len(qs) > 15 else ""))
        else:
            print("no imaginary modes below -%.1f cm^-1 (min %.2f cm^-1)" % (args.imag_tol, w.min()))
        print("max frequency %.2f cm^-1 = %.3f THz = %.2f meV"
              % (w.max(), w.max() * CONV["THz"], w.max() * CONV["meV"]))
        if args.csv:
            np.savetxt(args.csv, np.column_stack([x, w * f]), delimiter=",", fmt="%.6f",
                       header="x," + ",".join("w%d_%s" % (i + 1, args.unit)
                                              for i in range(w.shape[1])), comments="")
            print("CSV written to %s" % args.csv)
    dos = None
    if args.dos:
        dos = np.loadtxt(args.dos, comments="#")
        e = dos[:, 0]
        de = np.gradient(e)
        print("phonon DOS: %d points, integral %.3f (should be ~3*nat)"
              % (len(e), np.sum(dos[:, 1] * de)))

    if args.output.lower() == "none" and not args.show:
        return
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ylab = {"cm-1": "Frequency (cm$^{-1}$)", "THz": "Frequency (THz)", "meV": "Energy (meV)"}[args.unit]
    if x is not None and dos is not None:
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(6.4, 4.0), sharey=True,
                                      gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05})
    elif x is not None:
        fig, ax = plt.subplots(figsize=(5.0, 4.0))
        ax2 = None
    else:
        fig, ax2 = plt.subplots(figsize=(4.0, 3.5))
        ax = None
    if ax is not None:
        ax.plot(x, w * f, color="#1f4e79", lw=1.1)
        ax.axhline(0, color="0.4", lw=0.5)
        for t in ticks:
            ax.axvline(t, color="0.6", lw=0.6)
        ax.set_xlim(x[0], x[-1])
        ax.set_xticks(ticks)
        if labels and len(labels) == len(ticks):
            ax.set_xticklabels([r"$\Gamma$" if s.upper() in ("G", "GAMMA") else s for s in labels])
        else:
            if labels:
                print("WARNING: %d labels for %d ticks - labels ignored" % (len(labels), len(ticks)),
                      file=sys.stderr)
            ax.set_xticklabels([])
        ax.set_ylabel(ylab)
        if args.ylim:
            ax.set_ylim(*args.ylim)
    if ax2 is not None:
        if ax is not None:
            ax2.plot(dos[:, 1], dos[:, 0] * f, color="k", lw=1.0)
            if args.pdos:
                for j in range(2, dos.shape[1]):
                    ax2.plot(dos[:, j], dos[:, 0] * f, lw=0.8, label="atom %d" % (j - 1))
                ax2.legend(frameon=False, fontsize=7)
            ax2.set_xlabel("DOS")
            ax2.set_xlim(left=0)
        else:
            ax2.plot(dos[:, 0] * f, dos[:, 1], color="k", lw=1.0, label="total")
            if args.pdos:
                for j in range(2, dos.shape[1]):
                    ax2.plot(dos[:, 0] * f, dos[:, j], lw=0.8, label="atom %d" % (j - 1))
            ax2.legend(frameon=False, fontsize=7)
            ax2.set_xlabel(ylab)
            ax2.set_ylabel("DOS (states/cm$^{-1}$, as in fldos)")
    fig.tight_layout()
    if args.output.lower() != "none":
        fig.savefig(args.output, dpi=args.dpi)
        print("figure written to %s" % args.output)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
