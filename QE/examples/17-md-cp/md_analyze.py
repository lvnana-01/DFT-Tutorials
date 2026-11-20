#!/usr/bin/env python3
"""Summarise a pw.x MD run (calculation='md') or a cp.x .evp file.

pw.x mode : parses "temperature", "Ekin + Etot (const)" and "!    total energy"
            lines from the pw.x output and prints mean T, its fluctuation and the
            drift of the conserved quantity (useful for NVE checks).
evp mode  : reads prefix.evp written by cp.x (columns: nfi time(ps) ekinc Tcell
            Tion etot enthal econs econt ...) and reports the drift of econs and
            econt and the mean fictitious electronic kinetic energy ekinc.
Optional  : --xyz writes the pw.x trajectory to extended XYZ with ASE.

Only reads text files; it does not run any electronic-structure code.

Examples:
  python md_analyze.py pw si64.md.out --skip 200
  python md_analyze.py pw si64.nve.out --xyz si64_nve.extxyz
  python md_analyze.py evp tmp/h2o.evp
"""
import argparse
import re
import sys

RY_TO_EV = 13.605693122994
HA_TO_EV = 2.0 * RY_TO_EV


def stats(values):
    n = len(values)
    if n == 0:
        return float("nan"), float("nan")
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / n
    return mean, var ** 0.5


def linear_slope(x, y):
    n = len(x)
    if n < 2:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    if sxx == 0.0:
        return float("nan")
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / sxx


def parse_pw(path):
    temps, consts, times = [], [], []
    t_re = re.compile(r"^\s*temperature\s*=\s*([-\d.Ee+]+)\s*K")
    c_re = re.compile(r"Ekin \+ Etot \(const\)\s*=\s*([-\d.Ee+]+)\s*Ry")
    time_re = re.compile(r"^\s*time\s*=\s*([-\d.Ee+]+)\s*pico-seconds")
    with open(path, errors="replace") as fh:
        for line in fh:
            m = time_re.search(line)
            if m:
                times.append(float(m.group(1)))
                continue
            m = t_re.search(line)
            if m:
                temps.append(float(m.group(1)))
                continue
            m = c_re.search(line)
            if m:
                consts.append(float(m.group(1)))
    return times, temps, consts


def cmd_pw(args):
    times, temps, consts = parse_pw(args.output)
    n = min(len(temps), len(consts))
    if n == 0:
        sys.exit("no MD steps found (is this a pw.x calculation='md' output?)")
    temps, consts = temps[args.skip:n], consts[args.skip:n]
    steps = list(range(args.skip, n))
    if len(times) >= n:
        xs = times[args.skip:n]
        unit = "ps"
    else:
        xs = [float(s) for s in steps]
        unit = "step"
    t_mean, t_std = stats(temps)
    slope = linear_slope(xs, consts)
    print(f"MD steps parsed      : {n} (skipping first {args.skip})")
    print(f"T mean +- std        : {t_mean:.1f} +- {t_std:.1f} K")
    if args.natoms:
        drift = slope * RY_TO_EV * 1000.0 / args.natoms
        print(f"conserved-E drift    : {drift:.4f} meV/atom/{unit}")
    else:
        print(f"conserved-E drift    : {slope:.3e} Ry/{unit} (give --natoms for meV/atom)")
    span = (max(consts) - min(consts)) * RY_TO_EV * 1000.0 if consts else float("nan")
    print(f"conserved-E range    : {span:.3f} meV (whole cell)")
    if args.xyz:
        try:
            from ase.io import read, write
        except ImportError:
            sys.exit("ASE is required for --xyz (pip install ase)")
        frames = read(args.output, index=":", format="espresso-out")
        write(args.xyz, frames, format="extxyz")
        print(f"wrote {len(frames)} frames to {args.xyz}")


def cmd_evp(args):
    rows = []
    with open(args.evp) as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            parts = line.split()
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    rows = rows[args.skip:]
    if not rows:
        sys.exit("no data rows in evp file")
    time = [r[1] for r in rows]
    ekinc = [r[2] for r in rows]
    tion = [r[4] for r in rows]
    econs = [r[7] for r in rows]
    econt = [r[8] for r in rows]
    print(f"rows parsed          : {len(rows)}")
    print(f"Tion mean +- std     : {stats(tion)[0]:.1f} +- {stats(tion)[1]:.1f} K")
    print(f"ekinc mean (max)     : {stats(ekinc)[0]:.3e} ({max(ekinc):.3e}) Ha")
    print(f"econs drift          : {linear_slope(time, econs) * HA_TO_EV * 1000:.4f} meV/ps (whole cell)")
    print(f"econt drift          : {linear_slope(time, econt) * HA_TO_EV * 1000:.4f} meV/ps (whole cell)")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="mode", required=True)
    pw = sub.add_parser("pw", help="analyse a pw.x MD output")
    pw.add_argument("output")
    pw.add_argument("--skip", type=int, default=0, help="discard the first N steps (equilibration)")
    pw.add_argument("--natoms", type=int, default=0, help="number of atoms, to report drift per atom")
    pw.add_argument("--xyz", default="", help="write the trajectory to this extended-XYZ file (needs ASE)")
    pw.set_defaults(func=cmd_pw)
    ev = sub.add_parser("evp", help="analyse a cp.x prefix.evp file")
    ev.add_argument("evp")
    ev.add_argument("--skip", type=int, default=0)
    ev.set_defaults(func=cmd_evp)
    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
