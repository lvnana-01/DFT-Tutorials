#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qe_parse_out.py -- summarise pw.x output files (scf / nscf / relax / vc-relax / md).
                   解析 pw.x 输出：总能、费米能或 HOMO/LUMO、SCF 历史、力、应力、磁矩、
                   优化后的最终结构、JOB DONE，可导出 JSON / CSV。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Lines that are parsed (format strings of PW/src in QE 7.x):
  "!    total energy              =   -15.84... Ry"      converged SCF energy (one per ionic step)
  "     total energy              =  ..."  + "estimated scf accuracy    <  ..."   SCF history
  "     iteration #  1     ecut=    30.00 Ry     beta= 0.70"
  "     convergence has been achieved in   8 iterations" / "convergence NOT achieved after"
  "     the Fermi energy is     6.4993 ev"  / "the spin up/dw Fermi energies are"
  "     highest occupied, lowest unoccupied level (ev):     6.2 6.8" / "highest occupied level (ev):"
  "     atom    1 type  1   force =     0.00 0.00 0.00"  (Ry/au) and "Total force ="
  "          total   stress  (Ry/bohr**3)                   (kbar)     P=       -1.23"
  "     total magnetization       =     2.20 Bohr mag/cell" / "absolute magnetization"
  "     atom    1 (R=0.400)  charge= 13.97  magn=  2.31"     (site moments, nspin=2)
  "     bfgs converged in ..." / "End of BFGS Geometry Optimization" / "Final energy =" ("Final enthalpy")
  "Begin final coordinates" ... "End final coordinates"   (CELL_PARAMETERS + ATOMIC_POSITIONS)
  "JOB DONE."
Energies are reported in Ry and eV (1 Ry = 13.605693122994 eV); forces in Ry/bohr and eV/A.

Usage examples
--------------
    python qe_parse_out.py si.scf.out
    python qe_parse_out.py relax.out --scf          # also print SCF history per ionic step
    python qe_parse_out.py relax.out --final-structure final.in   # write final cell+positions
    python qe_parse_out.py *.out --csv summary.csv  # one row per file
    python qe_parse_out.py fe.scf.out --json fe.json

Exit status: 0 if every file has "JOB DONE." and no unconverged SCF, else 1 (handy in scripts).
Dependencies: python >= 3.9 (standard library only).
"""
import argparse
import csv
import json
import os
import re
import sys

RY_EV = 13.605693122994
BOHR_A = 0.529177210903
NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"


def fnum(s):
    try:
        return float(s.replace("D", "E").replace("d", "e"))
    except ValueError:
        return float("nan")  # e.g. "********" overflow fields


def parse(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    r = {"file": path, "job_done": False, "energies_Ry": [], "scf_steps": [],
         "scf_converged": [], "fermi_eV": None, "homo_eV": None, "lumo_eV": None,
         "total_force": [], "forces": None, "pressure_kbar": [], "stress_kbar": None,
         "total_mag": [], "abs_mag": [], "site_mag": None, "relax": None,
         "final_structure": None, "nat": None, "ecutwfc": None, "ecutrho": None,
         "nks": None, "version": None, "wall_time": None, "errors": []}
    cur_scf = []
    i = 0
    n = len(lines)
    while i < n:
        ln = lines[i]
        s = ln.strip()
        if s.startswith("Program PWSCF"):
            m = re.search(r"v\.(\S+)", s)
            if m:
                r["version"] = m.group(1)
        elif s.startswith("number of atoms/cell"):
            r["nat"] = int(s.split("=")[1])
        elif s.startswith("kinetic-energy cutoff"):
            r["ecutwfc"] = fnum(s.split("=")[1].split()[0])
        elif s.startswith("charge density cutoff"):
            r["ecutrho"] = fnum(s.split("=")[1].split()[0])
        elif s.startswith("number of k points="):
            r["nks"] = int(s.split("=")[1].split()[0])
        elif s.startswith("iteration #"):
            cur_scf.append({"iter": int(s.split("#")[1].split()[0])})
        elif s.startswith("total energy") and "=" in s and cur_scf:
            cur_scf[-1]["E_Ry"] = fnum(s.split("=")[1].split()[0])
        elif s.startswith("estimated scf accuracy") and cur_scf:
            cur_scf[-1]["acc_Ry"] = fnum(s.split("<")[1].split()[0])
        elif s.startswith("convergence has been achieved"):
            r["scf_converged"].append(True)
        elif s.startswith("convergence NOT achieved"):
            r["scf_converged"].append(False)
        elif ln.startswith("!"):
            m = re.search(r"=\s*(" + NUM + r")\s*Ry", ln)
            if m:
                r["energies_Ry"].append(fnum(m.group(1)))
                r["scf_steps"].append(cur_scf)
                cur_scf = []
        elif s.startswith("the Fermi energy is"):
            r["fermi_eV"] = fnum(s.split()[4])
        elif s.startswith("the spin up/dw Fermi energies are"):
            v = re.findall(NUM, s)
            r["fermi_eV"] = [fnum(v[0]), fnum(v[1])]
        elif s.startswith("highest occupied, lowest unoccupied level"):
            v = re.findall(NUM, s.split(":")[1])
            r["homo_eV"], r["lumo_eV"] = fnum(v[0]), fnum(v[1])
        elif s.startswith("highest occupied level"):
            r["homo_eV"] = fnum(re.findall(NUM, s.split(":")[1])[0])
        elif s.startswith("Forces acting on atoms"):
            f = []
            j = i + 1
            while j < n and not lines[j].strip():
                j += 1
            while j < n and lines[j].strip().startswith("atom") and "force =" in lines[j]:
                v = lines[j].split("=")[1].split()
                f.append([fnum(x) for x in v[:3]])
                j += 1
            r["forces"] = f
            i = j
            continue
        elif s.startswith("Total force ="):
            r["total_force"].append(fnum(s.split("=")[1].split()[0]))
        elif "total   stress" in s and "P=" in s:
            r["pressure_kbar"].append(fnum(s.split("P=")[1].split()[0]))
            st = []
            for k in range(1, 4):
                if i + k < n:
                    v = lines[i + k].split()
                    if len(v) >= 6:
                        st.append([fnum(x) for x in v[3:6]])
            r["stress_kbar"] = st
        elif s.startswith("total magnetization"):
            v = re.findall(NUM, s.split("=")[1])
            if v:  # collinear: one number; noncollinear: three components
                r["total_mag"].append(fnum(v[0]) if len(v) == 1 else [fnum(x) for x in v[:3]])
        elif s.startswith("absolute magnetization"):
            r["abs_mag"].append(fnum(s.split("=")[1].split()[0]))
        elif s.startswith("Magnetic moment per site"):
            sm = []
            j = i + 1
            while j < n and lines[j].strip().startswith("atom") and "magn=" in lines[j]:
                m = re.search(r"charge=\s*(" + NUM + r")\s+magn=\s*(" + NUM + ")", lines[j])
                if m:
                    sm.append({"charge": fnum(m.group(1)), "magn": fnum(m.group(2))})
                j += 1
            r["site_mag"] = sm
        elif s.startswith("bfgs converged in"):
            r["relax"] = {"converged": True, "line": s}
        elif s.startswith("bfgs failed after"):
            r["relax"] = {"converged": False, "line": s}
        elif s.startswith("The maximum number of steps has been reached"):
            r["relax"] = {"converged": False, "line": s}
        elif s.startswith("Final energy") or s.startswith("Final enthalpy"):
            if r["relax"] is None:
                r["relax"] = {"converged": True, "line": ""}
            r["relax"]["final_Ry"] = fnum(s.split("=")[1].split()[0])
        elif s == "Begin final coordinates":
            block = []
            j = i + 1
            while j < n and lines[j].strip() != "End final coordinates":
                block.append(lines[j])
                j += 1
            r["final_structure"] = block
            i = j
        elif s.startswith("PWSCF") and "WALL" in s:
            r["wall_time"] = s.split(":", 1)[1].strip() if ":" in s else s
        elif s.startswith("JOB DONE"):
            r["job_done"] = True
        elif s.startswith("%%%%%%"):
            msg = []
            j = i + 1
            while j < n and not lines[j].strip().startswith("%%%%%%"):
                if lines[j].strip():
                    msg.append(lines[j].strip())
                j += 1
            r["errors"].append(" | ".join(msg))
            i = j
        i += 1
    if cur_scf:
        r["scf_steps"].append(cur_scf)  # unfinished SCF (no "!" line)
    return r


def report(r, show_scf):
    p = print
    p("=" * 72)
    p("file      : %s   (QE %s)" % (r["file"], r["version"] or "?"))
    p("status    : %s" % ("JOB DONE" if r["job_done"] else "NOT finished (no 'JOB DONE.')"))
    if r["errors"]:
        for e in r["errors"]:
            p("ERROR     : %s" % e)
    p("nat %s | ecutwfc %s Ry | ecutrho %s Ry | k points %s"
      % (r["nat"], r["ecutwfc"], r["ecutrho"], r["nks"]))
    nconv = r["scf_converged"]
    if nconv:
        p("SCF       : %d cycles, %d converged, %d NOT converged"
          % (len(nconv), sum(nconv), len(nconv) - sum(nconv)))
    if r["energies_Ry"]:
        e = r["energies_Ry"][-1]
        p("E_total   : %.8f Ry = %.6f eV   (%d ionic step(s))" % (e, e * RY_EV, len(r["energies_Ry"])))
        if r["nat"]:
            p("E/atom    : %.6f eV" % (e * RY_EV / r["nat"]))
    if r["fermi_eV"] is not None:
        p("E_Fermi   : %s eV" % r["fermi_eV"])
    if r["homo_eV"] is not None:
        if r["lumo_eV"] is not None:
            p("HOMO/LUMO : %.4f / %.4f eV  -> gap %.4f eV (on this k mesh)"
              % (r["homo_eV"], r["lumo_eV"], r["lumo_eV"] - r["homo_eV"]))
        else:
            p("HOMO      : %.4f eV (add nbnd to get LUMO)" % r["homo_eV"])
    if r["total_force"]:
        p("Tot force : %.6f Ry/bohr (last)" % r["total_force"][-1])
    if r["forces"]:
        fmax = max((sum(c * c for c in f)) ** 0.5 for f in r["forces"])
        p("max |F|   : %.6f Ry/bohr = %.4f eV/A" % (fmax, fmax * RY_EV / BOHR_A))
    if r["pressure_kbar"]:
        p("Pressure  : %.2f kbar (last)" % r["pressure_kbar"][-1])
    if r["total_mag"]:
        p("Total mag : %s Bohr mag/cell; abs %s" % (r["total_mag"][-1],
                                                     r["abs_mag"][-1] if r["abs_mag"] else "?"))
    if r["site_mag"]:
        p("Site magn : " + " ".join("%.3f" % a["magn"] for a in r["site_mag"]))
    if r["relax"]:
        p("Relax     : %s %s" % ("converged" if r["relax"].get("converged") else "NOT converged",
                                 r["relax"].get("line", "")))
        if "final_Ry" in r["relax"]:
            p("Final E/H : %.10f Ry" % r["relax"]["final_Ry"])
    if r["wall_time"]:
        p("Time      : %s" % r["wall_time"])
    if show_scf:
        for k, steps in enumerate(r["scf_steps"], 1):
            p("-- SCF of ionic step %d --" % k)
            for st in steps:
                p("   iter %3d  E = %s Ry  acc = %s Ry"
                  % (st["iter"], st.get("E_Ry", "?"), st.get("acc_Ry", "?")))


def write_final(r, path):
    if not r["final_structure"]:
        sys.exit("no 'Begin final coordinates' block in %s" % r["file"])
    keep = [ln for ln in r["final_structure"]
            if not ln.strip().startswith(("new unit-cell volume", "density ="))]
    with open(path, "w", encoding="ascii", errors="replace", newline="\n") as fh:
        fh.write("! final structure from %s (paste into a new pw.x input)\n" % r["file"])
        fh.write("\n".join(keep).strip() + "\n")
    print("final structure written to %s" % path)
    if any("alat=" in ln for ln in keep):
        print("NOTE: CELL_PARAMETERS (alat=...) - keep celldm(1)/the same alat in the new input.")


def main():
    ap = argparse.ArgumentParser(description="Summarise pw.x output files.")
    ap.add_argument("outputs", nargs="+", help="pw.x output file(s)")
    ap.add_argument("--scf", action="store_true", help="print SCF history")
    ap.add_argument("--json", help="write full parsed data to JSON")
    ap.add_argument("--csv", help="write one summary row per file to CSV")
    ap.add_argument("--final-structure", metavar="FILE",
                    help="write final CELL_PARAMETERS/ATOMIC_POSITIONS (relax) to FILE")
    ap.add_argument("-q", "--quiet", action="store_true", help="no text report")
    args = ap.parse_args()

    results = []
    for f in args.outputs:
        if not os.path.isfile(f):
            print("skip %s: not found" % f, file=sys.stderr)
            continue
        r = parse(f)
        results.append(r)
        if not args.quiet:
            report(r, args.scf)
    if not results:
        sys.exit(1)
    if args.final_structure:
        write_final(results[-1], args.final_structure)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(results if len(results) > 1 else results[0], fh, indent=1)
        print("JSON written to %s" % args.json)
    if args.csv:
        cols = ["file", "job_done", "nat", "ecutwfc", "ecutrho", "nks", "E_Ry", "E_eV",
                "E_per_atom_eV", "fermi_eV", "homo_eV", "lumo_eV", "pressure_kbar",
                "total_force", "total_mag", "n_ionic", "scf_not_converged"]
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for r in results:
                e = r["energies_Ry"][-1] if r["energies_Ry"] else None
                w.writerow([r["file"], r["job_done"], r["nat"], r["ecutwfc"], r["ecutrho"], r["nks"],
                            e, e * RY_EV if e is not None else None,
                            e * RY_EV / r["nat"] if (e is not None and r["nat"]) else None,
                            r["fermi_eV"], r["homo_eV"], r["lumo_eV"],
                            r["pressure_kbar"][-1] if r["pressure_kbar"] else None,
                            r["total_force"][-1] if r["total_force"] else None,
                            r["total_mag"][-1] if r["total_mag"] else None,
                            len(r["energies_Ry"]), r["scf_converged"].count(False)])
        print("CSV written to %s" % args.csv)
    bad = [r for r in results if not r["job_done"] or False in r["scf_converged"]]
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
