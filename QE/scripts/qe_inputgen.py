#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qe_inputgen.py -- write a pw.x input (ibrav=0, explicit units) from a structure file.
                  从 CIF/POSCAR/xyz 等结构生成 pw.x 输入：ibrav=0 + CELL_PARAMETERS angstrom、
                  按 SSSP JSON 设置赝势文件名与 ecutwfc/ecutrho、按 k 间距生成 K_POINTS automatic、
                  可选自旋极化（nspin=2 + starting_magnetization）与 crystal_b 能带路径。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

What it writes
--------------
  &CONTROL  calculation, prefix, outdir, pseudo_dir (+ tprnfor/tstress, forc_conv_thr ...)
  &SYSTEM   ibrav = 0, nat, ntyp, ecutwfc, ecutrho, occupations/smearing/degauss,
            nspin/starting_magnetization, nbnd (bands)
  &ELECTRONS conv_thr, mixing_beta
  &IONS / &CELL  only for relax / vc-relax
  ATOMIC_SPECIES   label  mass  pseudo-file
  CELL_PARAMETERS angstrom
  ATOMIC_POSITIONS crystal
  K_POINTS automatic  n1 n2 n3 0 0 0     (scf/relax/vc-relax/nscf)
  K_POINTS crystal_b                     (--calculation bands, path from SeeK-path or ASE)

SSSP JSON (e.g. SSSP_1.3.0_PBE_efficiency.json from Materials Cloud) maps each element to
{"filename": ..., "cutoff_wfc": ..., "cutoff_rho": ...}.  With --sssp the script uses the
largest cutoff_wfc / cutoff_rho over all elements of the structure (the SSSP rule for
compounds).  Without --sssp you must give --ecutwfc (and optionally --dual or --ecutrho)
and the pseudopotential names default to "<El>.UPF" (edit them!).

k-point grid from a spacing (same convention as VASP KSPACING, 2*pi included):
    n_i = max(1, ceil(2*pi*|b_i| / kspacing)),   b_i = reciprocal vectors in 1/Angstrom
    (|b_i| without the 2*pi factor);  e.g. --kspacing 0.2  ~ dense, 0.3 ~ medium.

Usage examples
--------------
    python qe_inputgen.py Si.cif --sssp SSSP_1.3.0_PBE_efficiency.json --kspacing 0.2 -o si.scf.in
    python qe_inputgen.py POSCAR --sssp sssp.json --calculation vc-relax --smearing mv --degauss 0.02
    python qe_inputgen.py Fe.cif --sssp sssp.json --magmom Fe=0.5 --smearing mv --degauss 0.02
    python qe_inputgen.py Si.cif --sssp sssp.json --calculation bands --kpath seekpath --npoints 100
    python qe_inputgen.py MoS2.vasp --format vasp --ecutwfc 50 --dual 8 --kgrid 12 12 1 --pseudo-dir ./pseudo

Dependencies: python >= 3.9, ase (numpy).  Optional: seekpath (--kpath seekpath, uses spglib).
Always check the result with qe_check_input.py and by reading it before running pw.x.
"""
import argparse
import json
import math
import os
import sys

BOHR = 0.529177210903  # Angstrom


def need(module, hint):
    try:
        return __import__(module)
    except ImportError:
        sys.exit("qe_inputgen.py needs '%s': %s" % (module, hint))


def read_structure(path, fmt):
    need("ase", "pip install ase")
    from ase.io import read
    try:
        atoms = read(path, format=fmt) if fmt else read(path)
    except Exception as exc:  # noqa: BLE001
        sys.exit("cannot read %s with ASE: %s (try --format)" % (path, exc))
    if not atoms.cell.any() or abs(atoms.cell.volume) < 1e-6:
        sys.exit("structure has no cell; pw.x needs a periodic cell (add a box first)")
    return atoms


def species_order(symbols):
    seen = []
    for s in symbols:
        if s not in seen:
            seen.append(s)
    return seen


def load_sssp(path, elements):
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    missing = [e for e in elements if e not in data]
    if missing:
        sys.exit("elements %s not found in %s" % (", ".join(missing), path))
    pseudos = {e: data[e]["filename"] for e in elements}
    ecutwfc = max(float(data[e]["cutoff_wfc"]) for e in elements)
    ecutrho = max(float(data[e]["cutoff_rho"]) for e in elements)
    return pseudos, ecutwfc, ecutrho


def kgrid_from_spacing(atoms, spacing):
    import numpy as np
    rec = atoms.cell.reciprocal()  # 1/Angstrom, no 2*pi
    out = []
    for i in range(3):
        n = int(math.ceil(2.0 * math.pi * float(np.linalg.norm(rec[i])) / spacing))
        out.append(max(1, n))
    return out


def kpath_seekpath(atoms):
    """Return (atoms_in_seekpath_primitive_cell, [(label, kfrac), ...] segments)."""
    seekpath = need("seekpath", "pip install seekpath")
    from ase import Atoms
    cell = (atoms.cell.array.tolist(), atoms.get_scaled_positions().tolist(),
            atoms.get_atomic_numbers().tolist())
    res = seekpath.get_path(cell, with_time_reversal=True)
    prim = Atoms(numbers=res["primitive_types"], cell=res["primitive_lattice"],
                 scaled_positions=res["primitive_positions"], pbc=True)
    coords = res["point_coords"]
    pts = []
    for a, b in res["path"]:
        if pts and pts[-1][0] == a:
            pts.append((b, coords[b]))
        else:
            if pts:
                pts[-1] = (pts[-1][0], pts[-1][1], "break")
            pts.append((a, coords[a]))
            pts.append((b, coords[b]))
    return prim, pts


def kpath_ase(atoms):
    bp = atoms.cell.bandpath()
    special = bp.special_points
    pts = []
    for seg in bp.path.split(","):
        labels = []
        i = 0
        while i < len(seg):  # labels like G, X, K, U, W, L, S1, Y1 ...
            j = i + 1
            while j < len(seg) and seg[j].isdigit():
                j += 1
            labels.append(seg[i:j])
            i = j
        if pts:
            pts[-1] = (pts[-1][0], pts[-1][1], "break")
        for lab in labels:
            pts.append((lab, list(special[lab])))
    return atoms, pts


def fmt_label(lab):
    return "G" if lab in ("G", "GAMMA", "\\Gamma", "Γ") else lab


def build_input(args):
    atoms = read_structure(args.structure, args.format)
    kpts_block = []
    if args.calculation == "bands":
        if args.kpath == "seekpath":
            atoms, pts = kpath_seekpath(atoms)
            print("NOTE: SeeK-path standardised primitive cell is written; use the SAME cell "
                  "for the scf run (rerun this script with --calculation scf --kpath-cell).",
                  file=sys.stderr)
        else:
            atoms, pts = kpath_ase(atoms)
        kpts_block.append("K_POINTS crystal_b")
        kpts_block.append("%d" % len(pts))
        for p in pts:
            npt = 1 if (len(p) == 3) else args.npoints
            if p is pts[-1]:
                npt = 1
            k = p[1]
            kpts_block.append("  %12.8f %12.8f %12.8f %4d  ! %s"
                              % (k[0], k[1], k[2], npt, fmt_label(p[0])))
    elif args.kpath_cell and args.kpath == "seekpath":
        atoms, _ = kpath_seekpath(atoms)

    symbols = atoms.get_chemical_symbols()
    elements = species_order(symbols)
    if args.sssp:
        pseudos, ecw, ecr = load_sssp(args.sssp, elements)
    else:
        pseudos = {e: "%s.UPF" % e for e in elements}
        ecw, ecr = None, None
    if args.pseudo:
        for item in args.pseudo:
            el, _, name = item.partition("=")
            pseudos[el] = name
    if args.ecutwfc is not None:
        ecw = args.ecutwfc
    if ecw is None:
        sys.exit("give --sssp JSON or --ecutwfc")
    if args.ecutrho is not None:
        ecr = args.ecutrho
    elif args.dual is not None or ecr is None:
        ecr = ecw * (args.dual if args.dual is not None else 4.0)

    if args.calculation != "bands":
        if args.kgrid:
            grid = args.kgrid
        else:
            grid = kgrid_from_spacing(atoms, args.kspacing)
        kpts_block = ["K_POINTS automatic",
                      "  %d %d %d  %d %d %d" % (grid[0], grid[1], grid[2],
                                                 args.kshift[0], args.kshift[1], args.kshift[2])]

    from ase.data import atomic_masses, atomic_numbers
    magmom = {}
    for item in args.magmom or []:
        el, _, val = item.partition("=")
        magmom[el] = float(val)

    L = []
    calc = args.calculation
    L.append("! generated by qe_inputgen.py from %s -- check before running" % os.path.basename(args.structure))
    L.append("&CONTROL")
    L.append("  calculation = '%s'" % calc)
    L.append("  prefix      = '%s'" % args.prefix)
    L.append("  outdir      = '%s'" % args.outdir)
    L.append("  pseudo_dir  = '%s'" % args.pseudo_dir)
    if calc in ("scf", "relax", "vc-relax"):
        L.append("  tprnfor     = .true.")
        L.append("  tstress     = .true.")
    if calc in ("relax", "vc-relax"):
        L.append("  forc_conv_thr = %.1e" % args.forc_conv_thr)
        L.append("  etot_conv_thr = %.1e" % args.etot_conv_thr)
    L.append("/")
    L.append("&SYSTEM")
    L.append("  ibrav = 0")
    L.append("  nat   = %d" % len(atoms))
    L.append("  ntyp  = %d" % len(elements))
    L.append("  ecutwfc = %.1f" % ecw)
    L.append("  ecutrho = %.1f" % ecr)
    if args.smearing:
        L.append("  occupations = 'smearing'")
        L.append("  smearing    = '%s'" % args.smearing)
        L.append("  degauss     = %g" % args.degauss)
    elif calc == "nscf" and args.tetra:
        L.append("  occupations = 'tetrahedra_opt'")
    else:
        L.append("  occupations = 'fixed'")
    if magmom:
        L.append("  nspin = 2")
        for i, el in enumerate(elements, 1):
            if el in magmom:
                L.append("  starting_magnetization(%d) = %g" % (i, magmom[el]))
    if args.nbnd:
        L.append("  nbnd = %d" % args.nbnd)
    if args.vdw:
        L.append("  vdw_corr = '%s'" % args.vdw)
    L.append("/")
    L.append("&ELECTRONS")
    L.append("  conv_thr    = %.1e" % (args.conv_thr_per_atom * len(atoms)))
    L.append("  mixing_beta = %g" % args.mixing_beta)
    if calc == "bands" or calc == "nscf":
        L.append("  diago_full_acc = .true.")
    L.append("/")
    if calc in ("relax", "vc-relax"):
        L.append("&IONS")
        L.append("  ion_dynamics = 'bfgs'")
        L.append("/")
    if calc == "vc-relax":
        L.append("&CELL")
        L.append("  cell_dynamics = 'bfgs'")
        L.append("  press_conv_thr = %g" % args.press_conv_thr)
        if args.cell_dofree:
            L.append("  cell_dofree = '%s'" % args.cell_dofree)
        L.append("/")
    L.append("ATOMIC_SPECIES")
    for el in elements:
        L.append("  %-3s %10.4f  %s" % (el, atomic_masses[atomic_numbers[el]], pseudos[el]))
    L.append("CELL_PARAMETERS angstrom")
    for v in atoms.cell.array:
        L.append("  %16.10f %16.10f %16.10f" % tuple(v))
    L.append("ATOMIC_POSITIONS crystal")
    for s, p in zip(symbols, atoms.get_scaled_positions(wrap=True)):
        L.append("  %-3s %14.10f %14.10f %14.10f" % (s, p[0], p[1], p[2]))
    L.extend(kpts_block)
    text = "\n".join(L) + "\n"

    summary = ("elements %s | ecutwfc %.1f Ry ecutrho %.1f Ry (dual %.1f) | %s"
               % (" ".join(elements), ecw, ecr, ecr / ecw, kpts_block[1].strip()
                  if calc != "bands" else "%s path points" % kpts_block[1]))
    return text, summary, pseudos


def main():
    p = argparse.ArgumentParser(description="Generate a pw.x input (ibrav=0) from a structure file.",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("structure", help="structure file readable by ASE (cif, POSCAR, xyz with cell ...)")
    p.add_argument("--format", help="ASE format name if it cannot be guessed (e.g. vasp)")
    p.add_argument("-o", "--output", help="output file (default: <prefix>.<calculation>.in)")
    p.add_argument("--force", action="store_true", help="overwrite existing output")
    p.add_argument("--calculation", default="scf",
                   choices=["scf", "nscf", "bands", "relax", "vc-relax"])
    p.add_argument("--prefix", default="pwscf")
    p.add_argument("--outdir", default="./tmp")
    p.add_argument("--pseudo-dir", default="./pseudo")
    p.add_argument("--sssp", help="SSSP JSON file: pseudo file names + recommended cutoffs")
    p.add_argument("--pseudo", nargs="+", metavar="El=FILE", help="override pseudo file names")
    p.add_argument("--ecutwfc", type=float, help="override ecutwfc (Ry)")
    p.add_argument("--ecutrho", type=float, help="override ecutrho (Ry)")
    p.add_argument("--dual", type=float, help="ecutrho = dual * ecutwfc (4 for NC, 8-12 for US/PAW)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--kspacing", type=float, default=0.25,
                   help="k spacing in 1/Angstrom (2*pi included) for K_POINTS automatic")
    g.add_argument("--kgrid", type=int, nargs=3, metavar=("N1", "N2", "N3"), help="explicit grid")
    p.add_argument("--kshift", type=int, nargs=3, default=[0, 0, 0], metavar=("S1", "S2", "S3"))
    p.add_argument("--kpath", choices=["seekpath", "ase"], default="seekpath",
                   help="source of the crystal_b path for --calculation bands")
    p.add_argument("--kpath-cell", action="store_true",
                   help="for scf/nscf: write the SeeK-path standardised primitive cell "
                        "(so that a later bands run uses the same cell)")
    p.add_argument("--npoints", type=int, default=40, help="points per segment in crystal_b")
    p.add_argument("--nbnd", type=int, help="number of bands (bands/nscf: add empty bands)")
    p.add_argument("--smearing", choices=["gaussian", "mv", "m-p", "fd"],
                   help="use occupations='smearing' with this smearing (metals)")
    p.add_argument("--degauss", type=float, default=0.02, help="smearing width (Ry)")
    p.add_argument("--tetra", action="store_true", help="nscf: occupations='tetrahedra_opt' (DOS)")
    p.add_argument("--magmom", nargs="+", metavar="El=VAL",
                   help="nspin=2 with starting_magnetization (-1..1) per element")
    p.add_argument("--vdw", choices=["grimme-d2", "grimme-d3", "dft-d3", "ts", "xdm", "mbd"],
                   help="vdw_corr value")
    p.add_argument("--conv-thr-per-atom", type=float, default=1e-9,
                   help="conv_thr = value * nat (Ry)")
    p.add_argument("--mixing-beta", type=float, default=0.4)
    p.add_argument("--forc-conv-thr", type=float, default=1e-4, help="Ry/bohr")
    p.add_argument("--etot-conv-thr", type=float, default=1e-5, help="Ry")
    p.add_argument("--press-conv-thr", type=float, default=0.5, help="kbar")
    p.add_argument("--cell-dofree", help="e.g. 'ibrav', 'shape', '2Dxy'")
    args = p.parse_args()

    text, summary, pseudos = build_input(args)
    out = args.output or "%s.%s.in" % (args.prefix, args.calculation)
    if os.path.exists(out) and not args.force:
        sys.exit("%s exists (use --force)" % out)
    with open(out, "w", encoding="ascii", newline="\n") as fh:
        fh.write(text)
    print("wrote %s: %s" % (out, summary))
    miss = [f for f in pseudos.values() if not os.path.isfile(os.path.join(args.pseudo_dir, f))]
    if miss:
        print("WARNING: not found in %s: %s" % (args.pseudo_dir, ", ".join(miss)), file=sys.stderr)


if __name__ == "__main__":
    main()
