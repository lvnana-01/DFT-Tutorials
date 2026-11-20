#!/usr/bin/env python3
"""Read key results from a pw.x data-file-schema.xml with the standard library only.

The XML written by pw.x (outdir/prefix.save/data-file-schema.xml, also copied to
outdir/prefix.xml) follows the official qes schema (QEF/qeschema). Energies are
in Hartree, lengths in Bohr, forces in Ha/Bohr, stress in Ha/Bohr^3.

Example:
  python parse_xml.py tmp/si.save/data-file-schema.xml
  python parse_xml.py tmp/si.xml --json result.json
"""
import argparse
import json
import sys
import xml.etree.ElementTree as ET

HA_TO_EV = 27.211386245988
BOHR_TO_ANG = 0.529177210903
HA_BOHR3_TO_KBAR = 2.9421015697e5  # 1 Ha/Bohr^3 = 29421.0 GPa = 2.9421e5 kbar


def floats(text):
    return [float(x) for x in text.split()]


def find(root, path):
    node = root.find(path)
    return None if node is None else node.text


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("xml", help="data-file-schema.xml or prefix.xml")
    p.add_argument("--json", default="", help="also write the summary to this JSON file")
    args = p.parse_args()

    root = ET.parse(args.xml).getroot()
    out = root.find("output")
    if out is None:
        sys.exit("no <output> element: the run did not finish or this is not a pw.x XML")

    res = {}
    conv = find(out, "convergence_info/scf_conv/convergence_achieved")
    res["scf_converged"] = conv is not None and conv.strip().lower() == "true"
    nsteps = find(out, "convergence_info/scf_conv/n_scf_steps")
    res["n_scf_steps"] = int(nsteps) if nsteps else None

    etot = find(out, "total_energy/etot")
    res["etot_Ha"] = float(etot)
    res["etot_eV"] = float(etot) * HA_TO_EV

    bs = out.find("band_structure")
    for tag in ("fermi_energy", "highestOccupiedLevel", "lowestUnoccupiedLevel"):
        val = bs.find(tag) if bs is not None else None
        if val is not None:
            res[tag + "_eV"] = float(val.text) * HA_TO_EV

    struct = out.find("atomic_structure")
    res["nat"] = int(struct.get("nat"))
    res["species"] = [a.get("name") for a in struct.findall("atomic_positions/atom")]
    res["positions_ang"] = [[x * BOHR_TO_ANG for x in floats(a.text)]
                            for a in struct.findall("atomic_positions/atom")]
    res["cell_ang"] = [[x * BOHR_TO_ANG for x in floats(struct.find("cell/" + v).text)]
                       for v in ("a1", "a2", "a3")]

    forces = out.find("forces")
    if forces is not None:
        f = floats(forces.text)
        res["forces_eV_per_ang"] = [[x * HA_TO_EV / BOHR_TO_ANG for x in f[i:i + 3]]
                                    for i in range(0, len(f), 3)]
    stress = out.find("stress")
    if stress is not None:
        s = floats(stress.text)
        res["stress_kbar"] = [[x * HA_BOHR3_TO_KBAR for x in s[i:i + 3]] for i in range(0, 9, 3)]
        # same convention as the "P=" line of pw.x: P = trace(sigma)/3
        res["pressure_kbar"] = (s[0] + s[4] + s[8]) / 3.0 * HA_BOHR3_TO_KBAR

    for key in ("scf_converged", "n_scf_steps", "etot_eV", "fermi_energy_eV",
                "highestOccupiedLevel_eV", "lowestUnoccupiedLevel_eV", "pressure_kbar"):
        if key in res:
            print(f"{key:26s}: {res[key]}")
    if "forces_eV_per_ang" in res:
        fmax = max(max(abs(c) for c in row) for row in res["forces_eV_per_ang"])
        print(f"{'max |F_i| (eV/A)':26s}: {fmax:.4e}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(res, fh, indent=2)
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
