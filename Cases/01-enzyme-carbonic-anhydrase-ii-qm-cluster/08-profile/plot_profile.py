#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""plot_profile.py -- build a relative free-energy table and a reaction profile
from a small CSV that YOU fill with your own computed numbers.

Part of DFT-Tutorials, Cases/01-enzyme-carbonic-anhydrase-ii-qm-cluster.

The shipped energies_template.csv contains only placeholders (NaN,
placeholder=yes).  This script refuses to plot placeholder rows, so it can
never draw invented numbers.

    G(state) = E_sp(def2-TZVP) + G_corr(def2-SVP)          [Hartree]
    dG       = (G(state) - G(reference)) * 627.5095        [kcal/mol]

Usage
-----
    python plot_profile.py energies.csv                  # table only
    python plot_profile.py energies.csv -o profile.png   # table + figure
    python plot_profile.py energies.csv --ref R --unit kj

Needs matplotlib only for the figure.
"""

import argparse
import csv
import math
import sys

HARTREE = {"kcal": 627.5095, "kj": 2625.4996}
UNIT_LABEL = {"kcal": "kcal/mol", "kj": "kJ/mol"}


def read_rows(path):
    with open(path, newline="") as fh:
        lines = [l for l in fh if l.strip() and not l.lstrip().startswith("#")]
    rows = list(csv.DictReader(lines))
    need = {"label", "state_type", "E_sp_hartree", "G_corr_hartree", "n_imag",
            "placeholder"}
    if not rows or not need.issubset(rows[0].keys()):
        sys.exit("ERROR: CSV must have columns %s" % ", ".join(sorted(need)))
    return rows


def check(rows):
    bad = []
    for r in rows:
        if r["placeholder"].strip().lower() != "no":
            bad.append("%s: still marked placeholder" % r["label"])
            continue
        for k in ("E_sp_hartree", "G_corr_hartree", "n_imag"):
            try:
                v = float(r[k])
            except ValueError:
                v = float("nan")
            if math.isnan(v):
                bad.append("%s: %s missing" % (r["label"], k))
        try:
            nimag = int(float(r["n_imag"]))
        except ValueError:
            continue
        want = 1 if r["state_type"].strip().lower() == "ts" else 0
        if nimag != want:
            bad.append("%s: %d imaginary frequencies, expected %d"
                       % (r["label"], nimag, want))
    return bad


def main(argv=None):
    ap = argparse.ArgumentParser(description="Relative free-energy profile from your CSV")
    ap.add_argument("csv")
    ap.add_argument("--ref", default=None, help="reference label (default: first row)")
    ap.add_argument("--unit", choices=("kcal", "kj"), default="kcal")
    ap.add_argument("-o", "--output", default=None, help="figure file (png/pdf/svg)")
    a = ap.parse_args(argv)

    rows = read_rows(a.csv)
    problems = check(rows)
    if problems:
        sys.stderr.write("Refusing to continue - the CSV is not filled with your results:\n")
        for p in problems:
            sys.stderr.write("  - %s\n" % p)
        return 2

    f = HARTREE[a.unit]
    G = dict((r["label"], float(r["E_sp_hartree"]) + float(r["G_corr_hartree"]))
             for r in rows)
    ref = a.ref or rows[0]["label"]
    if ref not in G:
        sys.exit("ERROR: reference %s not in CSV" % ref)

    print("%-8s %-8s %16s %14s %12s" % ("label", "type", "G (Hartree)",
                                         "dG (%s)" % UNIT_LABEL[a.unit], "n_imag"))
    for r in rows:
        lab = r["label"]
        print("%-8s %-8s %16.6f %14.1f %12s" % (lab, r["state_type"], G[lab],
                                                 (G[lab] - G[ref]) * f, r["n_imag"]))

    if a.output:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(5.0, 3.4))
        xs, ys = [], []
        for i, r in enumerate(rows):
            y = (G[r["label"]] - G[ref]) * f
            x0 = 2 * i
            ax.plot([x0 - 0.35, x0 + 0.35], [y, y], color="#1f4e79", lw=2.0,
                    solid_capstyle="round")
            ax.text(x0, y + 0.6, "%s\n%.1f" % (r["label"], y), ha="center",
                    va="bottom", fontsize=8, color="#222222")
            xs.append(x0)
            ys.append(y)
        for i in range(len(xs) - 1):
            ax.plot([xs[i] + 0.35, xs[i + 1] - 0.35], [ys[i], ys[i + 1]],
                    color="#8a8a8a", lw=1.0, ls="--")
        ax.set_xticks([])
        ax.set_ylabel("relative G (%s)" % UNIT_LABEL[a.unit])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.margins(y=0.25)
        fig.tight_layout()
        fig.savefig(a.output, dpi=300)
        print("figure written to %s" % a.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
