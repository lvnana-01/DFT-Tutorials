#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""cut_cluster.py -- cut a QM cluster model of the hCA II active site from PDB 2CBA.

Part of DFT-Tutorials, Cases/01-enzyme-carbonic-anhydrase-ii-qm-cluster.

Reads PDB entry 2CBA (human carbonic anhydrase II, 1.54 A; Hakansson et al.,
J. Mol. Biol. 227, 1192 (1992), DOI 10.1016/0022-2836(92)90531-N) and builds
a small quantum-chemical cluster model of the active site:

  * Zn(II) ion;
  * the Zn-bound water (identified from the LINK records) -> hydroxide;
  * His94, His96, His119 cut between C-alpha and C-beta -> 4-methylimidazole
    (CB kept at its crystal position as the methyl carbon, CA replaced by H);
  * Thr199 cut between C-alpha and C-beta -> ethanol (CB, OG1, CG2);
  * optionally Glu106 cut between CB and CG -> acetate (CG, CD, OE1, OE2);
  * optionally N second-shell waters closest to the Zn-OH- oxygen (only
    those within --water-cutoff, default 3.4 A; waters overlapping CO2 are
    dropped);
  * optionally CO2, either copied from PDB 3D92 (CO2-bound hCA II, 1.10 A;
    Domsic et al., J. Biol. Chem. 283, 30766 (2008), DOI
    10.1074/jbc.M805353200) after superposing the zinc ions, or placed
    geometrically.

Atom order (1-based, fixed so that input templates can refer to it):
    1 Zn, 2 O(hydroxide), 3 H(hydroxide), then His94, His96, His119, Thr199,
    [Glu106], [waters], and CO2 LAST (C, O1, O2).

Sizes: default (1 water, no Glu106) 51 atoms; + CO2 54 atoms;
       --with-glu106 adds 7 atoms (58 / 61 atoms).

Every hydrogen is PLACED BY A GEOMETRIC RULE.  These are starting guesses
only; the structure must be pre-optimised (xTB) and optimised (DFT) before
any geometry or energy is reported.  The script runs no quantum-chemistry
program; it only parses PDB text and writes text.  Standard library only.

Usage
-----
    python cut_cluster.py --pdb 2CBA.pdb --outdir build
    python cut_cluster.py --pdb 2CBA.pdb --outdir build_co2 --add-co2 --co2-pdb 3D92.pdb
    python cut_cluster.py --pdb 2CBA.pdb --outdir build_glu --with-glu106 --add-co2 --co2-pdb 3D92.pdb
    python cut_cluster.py --download --outdir build      # fetches 2CBA itself

Outputs (in --outdir)
---------------------
    cluster.xyz                   XYZ geometry (also used by xTB and ORCA)
    cluster_coords.txt            "charge mult" line + Gaussian coordinate block
    cluster_freeze.txt            Gaussian Opt=ModRedundant lines "X n F"
    cluster_constraints_orca.txt  ORCA %geom Constraints lines "{ C n C }" (0-based)
    cluster_fix_xtb.txt           xTB "$fix atoms: ..." block (1-based)
    atoms_index.csv               index, residue, atom name, element, note
    build_report.txt              what was kept/cut/added, charge, electrons,
                                  reaction-coordinate atom numbers
"""

import argparse
import math
import os
import sys
import urllib.request

# --------------------------------------------------------------------------
# vector helpers
# --------------------------------------------------------------------------


def vsub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vadd(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vscale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    n = norm(a)
    if n < 1e-8:
        raise ValueError("zero-length vector")
    return vscale(a, 1.0 / n)


def dist(a, b):
    return norm(vsub(a, b))


def perp_towards(w, ref):
    """Unit vector perpendicular to w, pointing as much as possible to ref."""
    w = unit(w)
    r = vsub(ref, vscale(w, dot(ref, w)))
    if norm(r) < 1e-6:
        trial = (1.0, 0.0, 0.0) if abs(w[0]) < 0.9 else (0.0, 1.0, 0.0)
        r = cross(w, trial)
    return unit(r)


# --------------------------------------------------------------------------
# PDB parsing
# --------------------------------------------------------------------------


class Atom(object):
    __slots__ = ("name", "altloc", "resname", "chain", "resseq", "xyz")

    def __init__(self, name, altloc, resname, chain, resseq, xyz):
        self.name = name
        self.altloc = altloc
        self.resname = resname
        self.chain = chain
        self.resseq = resseq
        self.xyz = xyz


def parse_pdb(path):
    """Return (atoms, links); a link is ((res, chain, seq, atom), (...))."""
    atoms = []
    links = []
    with open(path, "r") as fh:
        for line in fh:
            rec = line[:6]
            if rec in ("ATOM  ", "HETATM"):
                atoms.append(Atom(line[12:16].strip(), line[16],
                                  line[17:20].strip(), line[21],
                                  int(line[22:26]),
                                  (float(line[30:38]), float(line[38:46]),
                                   float(line[46:54]))))
            elif rec == "LINK  ":
                links.append(((line[17:20].strip(), line[21], int(line[22:26]),
                               line[12:16].strip()),
                              (line[47:50].strip(), line[51], int(line[52:56]),
                               line[42:46].strip())))
    return atoms, links


def pick(atoms, resname, resseq, name, chain):
    """First matching atom; blank altloc preferred, then 'A'."""
    cand = [a for a in atoms if a.resname == resname and a.resseq == resseq
            and a.name == name and a.chain == chain]
    cand.sort(key=lambda a: (a.altloc != " ", a.altloc))
    if not cand:
        raise KeyError("atom not found: %s %d %s chain %s"
                       % (resname, resseq, name, chain))
    return cand[0].xyz


# --------------------------------------------------------------------------
# hydrogen placement rules (starting guesses only)
# --------------------------------------------------------------------------


def h_methyl(c, nbr, ref, ch=1.09):
    """Three tetrahedral H on carbon c bonded to nbr; the first H is rotated
    towards ref (use the cut atom, e.g. CA, so that it plays the link-H role)."""
    w = unit(vsub(c, nbr))
    p1 = perp_towards(w, vsub(ref, c))
    p2 = cross(w, p1)
    alpha = math.radians(180.0 - 109.47)
    out = []
    for i in range(3):
        t = math.radians(120.0 * i)
        d = vadd(vscale(w, math.cos(alpha)),
                 vadd(vscale(p1, math.sin(alpha) * math.cos(t)),
                      vscale(p2, math.sin(alpha) * math.sin(t))))
        out.append(vadd(c, vscale(unit(d), ch)))
    return out


def h_methylene(c, n1, n2, ch=1.09):
    """Two tetrahedral H on carbon c bonded to n1 and n2."""
    u1 = unit(vsub(n1, c))
    u2 = unit(vsub(n2, c))
    bis = unit(vscale(vadd(u1, u2), -1.0))
    n = unit(cross(u1, u2))
    half = math.radians(54.75)
    return [vadd(c, vscale(unit(vadd(vscale(bis, math.cos(half)),
                                     vscale(n, s * math.sin(half)))), ch))
            for s in (+1.0, -1.0)]


def h_ring(x, n1, n2, length):
    """Ring C-H or N-H: in plane, opposite to the bisector of the ring bonds."""
    d = unit(vscale(vadd(unit(vsub(n1, x)), unit(vsub(n2, x))), -1.0))
    return vadd(x, vscale(d, length))


def h_pointing(o, target, oh=0.97):
    """O-H bond pointing from o towards target (H-bond donor guess)."""
    return vadd(o, vscale(unit(vsub(target, o)), oh))


def water_h(o, acceptor, oh=0.96, angle=104.5):
    """Water hydrogens: one H points at the acceptor, the second one closes
    the H-O-H angle (orientation is arbitrary; xTB will fix it)."""
    a = unit(vsub(acceptor, o))
    p = perp_towards(a, (0.0, 0.0, 1.0))
    ang = math.radians(angle)
    h1 = vadd(o, vscale(a, oh))
    h2 = vadd(o, vscale(unit(vadd(vscale(a, math.cos(ang)),
                                  vscale(p, math.sin(ang)))), oh))
    return [h1, h2]


# --------------------------------------------------------------------------
# model construction
# --------------------------------------------------------------------------

HIS_RING = ("CG", "ND1", "CD2", "CE1", "NE2")
HIS_NBRS = {"CG": ("ND1", "CD2"), "ND1": ("CG", "CE1"), "CD2": ("CG", "NE2"),
            "CE1": ("ND1", "NE2"), "NE2": ("CD2", "CE1")}


class Model(object):
    def __init__(self):
        self.atoms = []     # list of dicts
        self.frozen = []    # 0-based indices

    def add(self, el, xyz, res, seq, name, note="", frozen=False):
        self.atoms.append(dict(el=el, xyz=xyz, res=res, seq=seq, name=name,
                               note=note))
        if frozen:
            self.frozen.append(len(self.atoms) - 1)
        return len(self.atoms) - 1


def build(args, atoms, links):
    ch = args.chain
    m = Model()
    rep = []

    def g(res, seq, name):
        return pick(atoms, res, seq, name, ch)

    # Zn-coordinating His nitrogens and the Zn-bound water, from LINK records
    zn_n = {}
    wat = None
    for a, b in links:
        for x, y in ((a, b), (b, a)):
            if x[0] == "HIS" and y[0] == "ZN":
                zn_n[x[2]] = x[3]
            if x[0] == "HOH" and y[0] == "ZN":
                wat = x[2]
    for seq in (94, 96, 119):
        if seq not in zn_n:
            sys.exit("ERROR: no LINK record between HIS %d and ZN" % seq)
    if wat is None:
        sys.exit("ERROR: no LINK record between ZN and a water")
    rep.append("Zn ligands from LINK records: "
               + ", ".join("His%d %s" % (k, v) for k, v in sorted(zn_n.items()))
               + ", HOH %d O" % wat)

    zn = g("ZN", 262, "ZN")
    wo = g("HOH", wat, "O")
    og1 = g("THR", 199, "OG1")

    # atoms 1-3: Zn and hydroxide
    m.add("Zn", zn, "ZN", 262, "ZN", "zinc(II)")
    o_idx = m.add("O", wo, "OH", wat, "O", "Zn-bound hydroxide O")
    if args.oh_mode == "thr":
        hoh = h_pointing(wo, og1)
    else:
        hoh = vadd(wo, vscale(unit(vsub(wo, zn)), 0.97))
    m.add("H", hoh, "OH", wat, "H", "hydroxide H (guess)")
    rep.append("Zn-O(water) in crystal: %.2f A; HOH %d -> OH-" % (dist(zn, wo), wat))

    # histidines -> 4-methylimidazole
    for seq in (94, 96, 119):
        pos = dict((n, g("HIS", seq, n)) for n in HIS_RING)
        for n in HIS_RING:
            m.add("N" if n[0] == "N" else "C", pos[n], "HIS", seq, n, "ring")
        cb = g("HIS", seq, "CB")
        ca = g("HIS", seq, "CA")
        m.add("C", cb, "HIS", seq, "CB", "methyl C (cut CA-CB); frozen", frozen=True)
        for j, h in enumerate(h_methyl(cb, pos["CG"], ca)):
            m.add("H", h, "HIS", seq, "HB%d" % (j + 1),
                  "link H replacing CA" if j == 0 else "methyl H")
        for n in ("CD2", "CE1"):
            a, b = HIS_NBRS[n]
            m.add("H", h_ring(pos[n], pos[a], pos[b], 1.08), "HIS", seq, "H" + n[1:],
                  "ring C-H")
        prot = "ND1" if zn_n[seq] == "NE2" else "NE2"
        a, b = HIS_NBRS[prot]
        m.add("H", h_ring(pos[prot], pos[a], pos[b], 1.01), "HIS", seq, "H" + prot[1:],
              "tautomer H on non-coordinating N")
        rep.append("His%d binds Zn via %s -> H placed on %s" % (seq, zn_n[seq], prot))

    # Thr199 -> ethanol (CB frozen; one CB-H plays the link-H role for CA)
    cb = g("THR", 199, "CB")
    cg2 = g("THR", 199, "CG2")
    m.add("C", cb, "THR", 199, "CB", "C (cut CA-CB); frozen", frozen=True)
    m.add("O", og1, "THR", 199, "OG1", "hydroxyl O")
    m.add("C", cg2, "THR", 199, "CG2", "methyl C")
    for j, h in enumerate(h_methylene(cb, og1, cg2)):
        m.add("H", h, "THR", 199, "HB%d" % (j + 1), "H on CB (one replaces CA)")
    for j, h in enumerate(h_methyl(cg2, cb, og1)):
        m.add("H", h, "THR", 199, "HG2%d" % (j + 1), "methyl H")
    # in the crystal Thr199 OG1 donates to Glu106 OE1 (about 2.6 A)
    m.add("H", h_pointing(og1, g("GLU", 106, "OE1")), "THR", 199, "HG1",
          "OG1-H towards Glu106 (guess)")

    # Glu106 -> acetate (optional)
    if args.with_glu106:
        cg = g("GLU", 106, "CG")
        cd = g("GLU", 106, "CD")
        m.add("C", cg, "GLU", 106, "CG", "methyl C (cut CB-CG); frozen", frozen=True)
        m.add("C", cd, "GLU", 106, "CD", "carboxylate C")
        m.add("O", g("GLU", 106, "OE1"), "GLU", 106, "OE1", "carboxylate O")
        m.add("O", g("GLU", 106, "OE2"), "GLU", 106, "OE2", "carboxylate O")
        for j, h in enumerate(h_methyl(cg, cd, g("GLU", 106, "CB"))):
            m.add("H", h, "GLU", 106, "HG%d" % (j + 1),
                  "link H replacing CB" if j == 0 else "methyl H")

    # CO2 coordinates (appended at the very end)
    co2 = []
    if args.add_co2 and args.co2_pdb:
        b_atoms, _ = parse_pdb(args.co2_pdb)
        shift = vsub(zn, pick(b_atoms, "ZN", 262, "ZN", ch))
        for n in ("C", "O1", "O2"):
            co2.append((n, vadd(pick(b_atoms, "CO2", 301, n, ch), shift)))
        rep.append("CO2 A 301 copied from %s, translated by (%.2f %.2f %.2f) A "
                   "to superpose Zn" % (args.co2_pdb, shift[0], shift[1], shift[2]))
    elif args.add_co2:
        # geometric guess: C 2.8 A from O(OH-), Zn-O-C about 105 deg, CO2 axis
        # perpendicular to the O...C line
        e1 = unit(vsub(wo, zn))
        p = perp_towards(e1, vsub(og1, wo))
        th = math.radians(105.0)
        c = vadd(wo, vscale(unit(vadd(vscale(e1, math.cos(th)),
                                      vscale(p, -math.sin(th)))), 2.80))
        axis = unit(cross(e1, p))
        co2 = [("C", c), ("O1", vadd(c, vscale(axis, 1.16))),
               ("O2", vadd(c, vscale(axis, -1.16)))]
        rep.append("CO2 placed geometrically (prefer --co2-pdb 3D92.pdb)")

    # second-shell waters (optional)
    if args.waters > 0:
        cand = [a for a in atoms if a.resname == "HOH" and a.name == "O"
                and a.chain == ch and a.resseq != wat and a.altloc in (" ", "A")]
        if co2:
            keep = []
            for a in cand:
                if min(dist(a.xyz, x[1]) for x in co2) < 2.2:
                    rep.append("HOH %d dropped: overlaps with CO2 (deep-water site)"
                               % a.resseq)
                else:
                    keep.append(a)
            cand = keep
        cand.sort(key=lambda a: dist(a.xyz, wo))
        for a in cand[:args.waters]:
            d = dist(a.xyz, wo)
            if d > args.water_cutoff:
                rep.append("no further water within %.2f A of OH-; stopping"
                           % args.water_cutoff)
                break
            m.add("O", a.xyz, "HOH", a.resseq, "O",
                  "second-shell water %.2f A from OH-; frozen" % d, frozen=True)
            for j, h in enumerate(water_h(a.xyz, wo)):
                m.add("H", h, "HOH", a.resseq, "H%d" % (j + 1), "water H (guess)")
            rep.append("second-shell water HOH %d kept (O...O(H) %.2f A)" % (a.resseq, d))

    c_idx = None
    for n, xyz in co2:
        i = m.add("C" if n == "C" else "O", xyz, "CO2", 301, n, "substrate")
        if n == "C":
            c_idx = i
    if c_idx is not None:
        rep.append("REACTION COORDINATE: O(OH-) = atom %d, C(CO2) = atom %d "
                   "(1-based); ORCA 0-based: %d and %d; start distance %.2f A"
                   % (o_idx + 1, c_idx + 1, o_idx, c_idx,
                      dist(m.atoms[o_idx]["xyz"], m.atoms[c_idx]["xyz"])))
    return m, rep


def charge_mult(args):
    """Zn(II) +2, OH- -1, 3 neutral His, neutral Thr199, neutral waters/CO2,
    Glu106 carboxylate -1 if present.  Zn(II) is d10 -> singlet."""
    q = 2 - 1
    if args.with_glu106:
        q -= 1
    return q, 1


# --------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------

ZNUM = {"H": 1, "C": 6, "N": 7, "O": 8, "Zn": 30}


def write(path, text):
    with open(path, "w", newline="\n") as fh:
        fh.write(text)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Cut a QM cluster model of the "
                                 "hCA II active site from PDB 2CBA.")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdb", help="local 2CBA.pdb")
    src.add_argument("--download", action="store_true",
                     help="download 2CBA.pdb from files.rcsb.org into --outdir")
    ap.add_argument("--outdir", default="build")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--with-glu106", action="store_true",
                    help="add Glu106 as acetate (total charge becomes 0)")
    ap.add_argument("--waters", type=int, default=1,
                    help="second-shell waters to keep (default 1; 0 = none)")
    ap.add_argument("--water-cutoff", type=float, default=3.4,
                    help="keep only waters whose O is within this distance (A) of "
                         "the hydroxide O (default 3.4)")
    ap.add_argument("--oh-mode", choices=("thr", "anti"), default="thr",
                    help="hydroxide H towards Thr199 OG1 (default) or along Zn-O")
    ap.add_argument("--add-co2", action="store_true", help="add CO2 (reactant complex)")
    ap.add_argument("--co2-pdb", default=None,
                    help="3D92.pdb; copy its CO2 A 301 after Zn superposition")
    args = ap.parse_args(argv)

    out = os.path.abspath(args.outdir)
    if not os.path.isdir(out):
        os.makedirs(out)
    if args.download:
        pdb = os.path.join(out, "2CBA.pdb")
        urllib.request.urlretrieve("https://files.rcsb.org/download/2CBA.pdb", pdb)
    else:
        pdb = args.pdb
        if not os.path.isfile(pdb):
            sys.exit("ERROR: cannot read %s" % pdb)

    atoms, links = parse_pdb(pdb)
    m, rep = build(args, atoms, links)
    q, mult = charge_mult(args)
    nel = sum(ZNUM[a["el"]] for a in m.atoms) - q

    title = ("hCA II cluster from 2CBA, charge %d mult %d, H placed geometrically, "
             "NOT optimised" % (q, mult))
    xyz = ["%d" % len(m.atoms), title]
    gau = ["%d %d" % (q, mult)]
    for a in m.atoms:
        x, y, z = a["xyz"]
        xyz.append("%-2s %12.6f %12.6f %12.6f" % (a["el"], x, y, z))
        gau.append("%-2s %14.6f %14.6f %14.6f" % (a["el"], x, y, z))
    write(os.path.join(out, "cluster.xyz"), "\n".join(xyz) + "\n")
    write(os.path.join(out, "cluster_coords.txt"), "\n".join(gau) + "\n")
    write(os.path.join(out, "cluster_freeze.txt"),
          "".join("X %d F\n" % (i + 1) for i in m.frozen))
    write(os.path.join(out, "cluster_constraints_orca.txt"),
          "".join("    { C %d C }\n" % i for i in m.frozen))
    write(os.path.join(out, "cluster_fix_xtb.txt"),
          "$fix\n   atoms: %s\n$end\n" % ",".join(str(i + 1) for i in m.frozen))
    rows = ["index,residue,resseq,atom,element,note"]
    for i, a in enumerate(m.atoms):
        rows.append("%d,%s,%d,%s,%s,%s" % (i + 1, a["res"], a["seq"], a["name"],
                                           a["el"], a["note"].replace(",", ";")))
    write(os.path.join(out, "atoms_index.csv"), "\n".join(rows) + "\n")

    lines = ["hCA II QM cluster build report",
             "source PDB      : %s" % pdb,
             "atoms           : %d" % len(m.atoms),
             "charge / mult   : %d / %d" % (q, mult),
             "electrons       : %d (%s)" % (nel, "even, closed shell OK"
                                            if nel % 2 == 0 else "ODD - check model"),
             "frozen (1-based): %s" % ", ".join(str(i + 1) for i in m.frozen),
             ""] + rep + [
             "",
             "charge bookkeeping: Zn(II) +2, OH- -1, His94/96/119 0, Thr199 0%s"
             % (", Glu106 -1" if args.with_glu106 else ""),
             "All hydrogens are geometric guesses. Pre-optimise (01-xtb-preopt) and",
             "optimise with DFT (02-opt) before reporting anything."]
    write(os.path.join(out, "build_report.txt"), "\n".join(lines) + "\n")
    sys.stdout.write("\n".join(lines) + "\n\nfiles written to %s\n" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
