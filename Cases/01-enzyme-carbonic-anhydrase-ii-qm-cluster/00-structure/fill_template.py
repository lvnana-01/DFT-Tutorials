#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fill_template.py -- paste coordinates / frozen-atom lines into an input template.

Part of DFT-Tutorials, Cases/01-enzyme-carbonic-anhydrase-ii-qm-cluster.

The Gaussian (.gjf) and ORCA (.inp) templates of this case study contain
placeholders instead of coordinates, because no calculation was run for the
tutorial and no optimised geometry is shipped.  This helper replaces

    @COORDS@        by the coordinate block of an XYZ file   (needs --xyz)
    @CHARGE_MULT@   by "charge multiplicity"                 (needs --charge --mult)
    @FREEZE@        by Gaussian ModRedundant lines "X n F"   (1-based)
    @ORCA_FIX@      by ORCA Constraints lines "{ C n C }"    (0-based)

Usage
-----
    python fill_template.py ../02-opt/g16_opt_reactant.gjf opt_reactant.gjf \
        --xyz ../01-xtb-preopt/xtb_reactant/xtbopt.xyz --charge 1 --mult 1 \
        --frozen-file build_co2/cluster_freeze.txt

    # templates that read the geometry from a .chk need no --xyz
    python fill_template.py ../03-ts/g16_scan_OC.gjf scan_OC.gjf \
        --frozen-file build_co2/cluster_freeze.txt

    python fill_template.py ../02-opt/orca_opt_reactant.inp orca_opt_reactant.inp \
        --frozen 9,21,33,40,49

The script stops if any @PLACEHOLDER@ is left.  Standard library only; no
calculation is performed.
"""

import argparse
import re
import sys


def read_xyz(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    body = []
    for line in lines[2:2 + n]:
        p = line.split()
        body.append("%-2s %14.6f %14.6f %14.6f"
                    % (p[0], float(p[1]), float(p[2]), float(p[3])))
    if len(body) != n:
        sys.exit("ERROR: %s announces %d atoms but has %d lines" % (path, n, len(body)))
    return body


def frozen_from_file(path):
    idx = []
    with open(path) as fh:
        for line in fh:
            m = re.match(r"\s*X\s+(\d+)\s+F", line)
            if m:
                idx.append(int(m.group(1)))
    return idx


def main(argv=None):
    ap = argparse.ArgumentParser(description="Fill @PLACEHOLDERS@ in an input template")
    ap.add_argument("template")
    ap.add_argument("output")
    ap.add_argument("--xyz", default=None, help="XYZ file used for @COORDS@")
    ap.add_argument("--charge", type=int, default=None)
    ap.add_argument("--mult", type=int, default=None)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--frozen", default="", help="comma-separated 1-based atom numbers")
    g.add_argument("--frozen-file", default=None,
                   help="cluster_freeze.txt written by cut_cluster.py")
    a = ap.parse_args(argv)

    if a.frozen_file:
        frozen = frozen_from_file(a.frozen_file)
    else:
        frozen = [int(x) for x in a.frozen.split(",") if x.strip()]

    with open(a.template) as fh:
        text = fh.read()

    rep = {"@FREEZE@": "\n".join("X %d F" % i for i in frozen),
           "@ORCA_FIX@": "\n".join("    { C %d C }" % (i - 1) for i in frozen)}
    natom = None
    if a.xyz:
        coords = read_xyz(a.xyz)
        natom = len(coords)
        for i in frozen:
            if i < 1 or i > natom:
                sys.exit("ERROR: frozen atom %d outside 1..%d" % (i, natom))
        rep["@COORDS@"] = "\n".join(coords)
    if a.charge is not None and a.mult is not None:
        rep["@CHARGE_MULT@"] = "%d %d" % (a.charge, a.mult)
    for k, v in rep.items():
        text = text.replace(k, v)

    left = re.findall(r"@[A-Z_]+@", text)
    if left:
        sys.exit("ERROR: unreplaced placeholders: %s" % ", ".join(sorted(set(left))))
    if a.output.endswith(".gjf"):
        text = text.rstrip("\n") + "\n\n"       # Gaussian needs a final blank line
    with open(a.output, "w", newline="\n") as fh:
        fh.write(text)
    sys.stdout.write("wrote %s (%s, %d frozen atoms)\n"
                     % (a.output, "%d atoms" % natom if natom else "no coordinates pasted",
                        len(frozen)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
