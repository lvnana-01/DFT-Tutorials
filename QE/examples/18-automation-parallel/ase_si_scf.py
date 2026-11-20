#!/usr/bin/env python3
"""Write (and optionally run) a pw.x SCF for bulk Si with ASE's Espresso calculator.

Uses the profile-based API of ASE >= 3.23 (EspressoProfile). The older
Espresso(command=...) form raises an error in recent ASE versions.

By default the script ONLY writes the input file (espresso.pwi) so that you can
inspect it; pass --run to actually launch pw.x through the profile command.

Examples:
  python ase_si_scf.py --pseudo-dir ./pseudo                 # write input only
  python ase_si_scf.py --pseudo-dir ./pseudo --run \
      --command "mpirun -np 4 pw.x -nk 2"                    # run pw.x
"""
import argparse
from pathlib import Path

from ase.build import bulk


def build_parameters():
    # Nested input_data dict: section names follow the pw.x namelists.
    input_data = {
        "control": {
            "calculation": "scf",
            "prefix": "si_ase",
            "outdir": "./tmp",
            "tprnfor": True,
            "tstress": True,
        },
        "system": {
            "ecutwfc": 30.0,
            "ecutrho": 240.0,
            "occupations": "fixed",
        },
        "electrons": {
            "conv_thr": 1.0e-10,
            "mixing_beta": 0.7,
        },
    }
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}  # SSSP 1.3 PBE efficiency
    return input_data, pseudopotentials


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pseudo-dir", default="./pseudo", help="directory with the UPF files")
    p.add_argument("--workdir", default="ase_si", help="calculation directory")
    p.add_argument("--command", default="pw.x", help="launcher + pw.x (+ parallel flags)")
    p.add_argument("--kpts", type=int, nargs=3, default=[8, 8, 8])
    p.add_argument("--run", action="store_true", help="run pw.x (default: write input only)")
    args = p.parse_args()

    atoms = bulk("Si", "diamond", a=5.43)  # Angstrom; use your converged lattice constant
    input_data, pseudos = build_parameters()

    if not args.run:
        from ase.io import write

        workdir = Path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        input_data["control"]["pseudo_dir"] = str(Path(args.pseudo_dir).resolve())
        out = workdir / "espresso.pwi"
        write(out, atoms, format="espresso-in", input_data=input_data,
              pseudopotentials=pseudos, kpts=tuple(args.kpts))
        print(f"wrote {out} (not executed; add --run to launch pw.x)")
        return

    from ase.calculators.espresso import Espresso, EspressoProfile

    profile = EspressoProfile(command=args.command, pseudo_dir=str(Path(args.pseudo_dir).resolve()))
    atoms.calc = Espresso(profile=profile, directory=args.workdir, input_data=input_data,
                          pseudopotentials=pseudos, kpts=tuple(args.kpts))
    energy = atoms.get_potential_energy()   # eV
    forces = atoms.get_forces()             # eV/Angstrom
    stress = atoms.get_stress()             # eV/Angstrom^3, Voigt order
    print(f"E = {energy:.6f} eV")
    print(f"max |F| = {abs(forces).max():.4e} eV/A")
    print(f"stress (Voigt, eV/A^3) = {stress}")


if __name__ == "__main__":
    main()
