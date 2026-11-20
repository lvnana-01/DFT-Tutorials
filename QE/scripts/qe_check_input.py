#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qe_check_input.py -- static sanity check of pw.x input files (no pw.x needed).
                     pw.x 输入静态检查：namelist 是否以 / 结束、CRLF/非 ASCII 字符、卡片单位缺失、
                     已删除/废弃关键字（lda_plus_u、Hubbard_U、wf_collect ...）、ecutrho/ecutwfc 比例、
                     nat/ntyp 与卡片行数是否一致、赝势文件是否存在、ibrav=0 是否给了 CELL_PARAMETERS。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Checks (ERROR = pw.x will stop, WARN = probably wrong or outdated, INFO = hint)
  ERROR  namelist (&CONTROL, &SYSTEM, &ELECTRONS, ...) not closed by "/"
  ERROR  keyword removed in QE 7.1 (old DFT+U syntax: lda_plus_u, lda_plus_u_kind, Hubbard_U(i),
         Hubbard_J0, Hubbard_V, Hubbard_J, U_projection_type)
         -> use the HUBBARD card, e.g. "HUBBARD {ortho-atomic}" / "U Fe-3d 4.0"
  ERROR  ibrav = 0 without CELL_PARAMETERS; both celldm(1) and A given
  ERROR  ATOMIC_SPECIES / ATOMIC_POSITIONS line count differs from ntyp / nat
  ERROR  ecutrho/ecutwfc < 4, K_POINTS automatic without 6 integers
  WARN   CRLF line endings (Windows editors), non-ASCII characters, TAB characters
  WARN   CELL_PARAMETERS / ATOMIC_POSITIONS without explicit unit (deprecated)
  WARN   obsolete keywords that are ignored: wf_collect, lkpoint_dir
  WARN   ecutrho = 4*ecutwfc with ultrasoft/PAW pseudopotentials (usually 8-12)
  WARN   no K_POINTS card
  WARN   occupations='smearing' without degauss, nspin=2 without starting_magnetization
  WARN   pseudopotential file not found in pseudo_dir ($ESPRESSO_PSEUDO, ~/espresso/pseudo/)
  INFO   calculation='bands' with K_POINTS automatic, relax without tprnfor, etc.

Usage examples
--------------
    python qe_check_input.py si.scf.in
    python qe_check_input.py *.in --pseudo-dir ./pseudo
    python qe_check_input.py fe.scf.in --strict      # exit 1 also on WARN

Exit status: 0 = no ERROR (and no WARN with --strict), 1 otherwise.
Dependencies: python >= 3.9 (standard library only).
"""
import argparse
import os
import re
import sys

NAMELISTS = ("CONTROL", "SYSTEM", "ELECTRONS", "IONS", "CELL", "FCP", "RISM")
CARDS = ("ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS",
         "CONSTRAINTS", "ATOMIC_VELOCITIES", "ATOMIC_FORCES", "ADDITIONAL_K_POINTS",
         "SOLVENTS", "HUBBARD")
REMOVED = {  # removed in QE 7.1 (new HUBBARD card)
    "lda_plus_u": "HUBBARD card", "lda_plus_u_kind": "HUBBARD card",
    "hubbard_u": "HUBBARD card: U <label>-<manifold> <value>",
    "hubbard_j0": "HUBBARD card: J0 ...", "hubbard_v": "HUBBARD card: V ...",
    "hubbard_j": "HUBBARD card: J ...", "u_projection_type": "HUBBARD {projection}",
}
OBSOLETE = {"wf_collect": "ignored since QE 6.x (always collected)",
            "lkpoint_dir": "ignored"}


class Checker:
    def __init__(self, path, pseudo_dir):
        self.path = path
        self.msgs = []
        self.pseudo_dir_cli = pseudo_dir

    def add(self, level, msg, line=None):
        self.msgs.append((level, msg, line))

    def run(self):
        raw = open(self.path, "rb").read()
        if b"\r\n" in raw:
            self.add("WARN", "CRLF line endings - convert with dos2unix (pw.x may misread values)")
        try:
            raw.decode("ascii")
        except UnicodeDecodeError:
            for k, ln in enumerate(raw.split(b"\n"), 1):
                try:
                    ln.decode("ascii")
                except UnicodeDecodeError:
                    self.add("WARN", "non-ASCII characters (full-width quotes/spaces?)", k)
                    break
        text = raw.decode("utf-8", errors="replace").replace("\r", "")
        lines = text.split("\n")
        if "\t" in text:
            self.add("INFO", "TAB characters present (usually harmless)")
        nl, kw = self.parse_namelists(lines)
        cards = self.parse_cards(lines)
        self.check_keywords(kw)
        self.check_cards(kw, cards, nl)
        return self.msgs

    def parse_namelists(self, lines):
        found = {}
        kw = {}
        cur = None
        for k, ln in enumerate(lines, 1):
            s = ln.split("!")[0].strip()
            if not s:
                continue
            if s.startswith("&"):
                if cur is not None:
                    self.add("ERROR", "namelist &%s not closed by '/' before &%s" % (cur, s[1:]), k)
                cur = s[1:].split()[0].upper()
                if cur not in NAMELISTS:
                    self.add("WARN", "unknown namelist &%s" % cur, k)
                found[cur] = k
                s = s[len(cur) + 1:].strip()
                if not s:
                    continue
            if cur is None:
                continue
            if s == "/" or s.endswith("/") and not re.search(r"['\"][^'\"]*/\s*$", s):
                body = s[:-1]
                self._kv(body, kw, k)
                cur = None
                continue
            if s.split()[0].upper() in CARDS:
                self.add("ERROR", "card %s found inside namelist &%s (missing '/')"
                         % (s.split()[0], cur), k)
                cur = None
                continue
            self._kv(s, kw, k)
        if cur is not None:
            self.add("ERROR", "namelist &%s not closed by '/'" % cur)
        for req in ("CONTROL", "SYSTEM", "ELECTRONS"):
            if req not in found:
                lvl = "WARN" if req == "CONTROL" else "ERROR"
                self.add(lvl, "namelist &%s missing (must be present even if empty)" % req)
        return found, kw

    def _kv(self, s, kw, k):
        for part in re.findall(r"([A-Za-z_][\w]*(?:\([\d,\s]+\))?)\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s]+)", s):
            key = part[0].lower().replace(" ", "")
            kw[key] = (part[1].strip("'\"").strip(), k)

    def parse_cards(self, lines):
        cards = {}
        cur = None
        for k, ln in enumerate(lines, 1):
            s = ln.split("!")[0].split("#")[0].strip()
            if not s:
                continue
            w = s.split()[0].upper()
            if w in CARDS:
                cur = w
                opt = s[len(w):].strip().strip("{}()").strip().lower()
                cards[cur] = {"line": k, "opt": opt, "rows": []}
                continue
            if s.startswith("&") or s == "/":
                cur = None
                continue
            if cur:
                cards[cur]["rows"].append((k, s))
        return cards

    @staticmethod
    def num(v):
        try:
            return float(v.lower().replace("d", "e"))
        except (ValueError, AttributeError):
            return None

    def check_keywords(self, kw):
        for key, (val, k) in kw.items():
            base = key.split("(")[0]
            if base in REMOVED:
                self.add("ERROR", "'%s' was removed in QE 7.1 -> use %s" % (key, REMOVED[base]), k)
            if base in OBSOLETE:
                self.add("WARN", "'%s' is obsolete: %s - delete it" % (key, OBSOLETE[base]), k)
        ew = self.num(kw.get("ecutwfc", (None,))[0])
        er = self.num(kw.get("ecutrho", (None,))[0])
        if ew is None:
            self.add("ERROR", "ecutwfc missing in &SYSTEM")
        elif er is not None:
            if er / ew < 4 - 1e-6:
                self.add("ERROR", "ecutrho/ecutwfc = %.2f < 4 (pw.x stops)" % (er / ew))
        occ = kw.get("occupations", ("",))[0].lower()
        if occ == "smearing" and "degauss" not in kw:
            self.add("WARN", "occupations='smearing' without degauss (default 0.D0 Ry)")
        nspin = kw.get("nspin", ("1",))[0]
        if nspin == "2" and not any(x.startswith("starting_magnetization") for x in kw) \
                and "tot_magnetization" not in kw:
            self.add("WARN", "nspin=2 without starting_magnetization/tot_magnetization "
                             "(pw.x stops or converges to non-magnetic state)")
        calc = kw.get("calculation", ("scf",))[0].lower()
        if calc in ("relax", "vc-relax", "md", "vc-md") and "tprnfor" not in kw:
            self.add("INFO", "forces are computed anyway for %s" % calc)
        if "celldm(1)" in kw and "a" in kw:
            self.add("ERROR", "both celldm(1) and A given - use only one")
        self.kw = kw
        self.ecut = (ew, er)

    def check_cards(self, kw, cards, nl):
        ibrav = kw.get("ibrav", (None,))[0]
        if ibrav is None:
            self.add("ERROR", "ibrav missing in &SYSTEM")
        elif ibrav.strip() == "0" and "CELL_PARAMETERS" not in cards:
            self.add("ERROR", "ibrav=0 but no CELL_PARAMETERS card")
        for c in ("CELL_PARAMETERS", "ATOMIC_POSITIONS"):
            if c in cards and not cards[c]["opt"]:
                self.add("WARN", "%s without unit (deprecated) - write e.g. %s %s"
                         % (c, c, "angstrom" if c == "CELL_PARAMETERS" else "crystal"),
                         cards[c]["line"])
        if "CELL_PARAMETERS" in cards:
            if len(cards["CELL_PARAMETERS"]["rows"]) < 3:
                self.add("ERROR", "CELL_PARAMETERS needs 3 rows", cards["CELL_PARAMETERS"]["line"])
            if cards["CELL_PARAMETERS"]["opt"] == "alat" and "celldm(1)" not in kw and "a" not in kw:
                self.add("ERROR", "CELL_PARAMETERS alat needs celldm(1) or A")
        nat = self.num(kw.get("nat", (None,))[0])
        ntyp = self.num(kw.get("ntyp", (None,))[0])
        species = []
        if "ATOMIC_SPECIES" not in cards:
            self.add("ERROR", "ATOMIC_SPECIES card missing")
        else:
            rows = cards["ATOMIC_SPECIES"]["rows"]
            species = [r[1].split() for r in rows]
            if ntyp is not None and len(rows) != int(ntyp):
                self.add("ERROR", "ntyp = %d but ATOMIC_SPECIES has %d lines" % (ntyp, len(rows)),
                         cards["ATOMIC_SPECIES"]["line"])
        if "ATOMIC_POSITIONS" not in cards:
            self.add("ERROR", "ATOMIC_POSITIONS card missing")
        else:
            rows = [r for r in cards["ATOMIC_POSITIONS"]["rows"]
                    if r[1].split()[0].upper() not in CARDS]
            if nat is not None and len(rows) != int(nat):
                self.add("ERROR", "nat = %d but ATOMIC_POSITIONS has %d lines" % (nat, len(rows)),
                         cards["ATOMIC_POSITIONS"]["line"])
            labels = {s[0] for s in species if s}
            for k, s in rows:
                if s.split()[0] not in labels and labels:
                    self.add("ERROR", "atom label '%s' not in ATOMIC_SPECIES" % s.split()[0], k)
                    break
        if "K_POINTS" not in cards:
            self.add("WARN", "no K_POINTS card - add e.g. K_POINTS automatic / K_POINTS gamma")
        else:
            opt = cards["K_POINTS"]["opt"]
            calc = kw.get("calculation", ("scf",))[0].lower()
            if calc == "bands" and opt in ("automatic", "gamma"):
                self.add("WARN", "calculation='bands' with K_POINTS %s - use crystal_b/tpiba_b" % opt)
            if opt in ("crystal_b", "tpiba_b") and calc not in ("bands", "nscf"):
                self.add("WARN", "K_POINTS %s in a %s run (meant for bands)" % (opt, calc))
            if opt == "automatic" and cards["K_POINTS"]["rows"]:
                v = cards["K_POINTS"]["rows"][0][1].split()
                if len(v) != 6:
                    self.add("ERROR", "K_POINTS automatic needs 'n1 n2 n3 s1 s2 s3'",
                             cards["K_POINTS"]["rows"][0][0])
        if "HUBBARD" in cards and not cards["HUBBARD"]["opt"]:
            self.add("WARN", "HUBBARD card without projection type, e.g. HUBBARD {ortho-atomic}")
        # pseudopotentials
        pdir = self.pseudo_dir_cli or kw.get("pseudo_dir", (None,))[0] \
            or os.environ.get("ESPRESSO_PSEUDO") or os.path.expanduser("~/espresso/pseudo")
        base = os.path.dirname(os.path.abspath(self.path))
        pdir_abs = pdir if os.path.isabs(pdir) else os.path.join(base, pdir)
        uspp = False
        for s in species:
            if len(s) < 3:
                self.add("ERROR", "ATOMIC_SPECIES line needs 'label mass file': %s" % " ".join(s))
                continue
            f = os.path.join(pdir_abs, s[2])
            if not os.path.isfile(f):
                self.add("WARN", "pseudopotential %s not found in %s" % (s[2], pdir_abs))
                continue
            try:
                head = open(f, encoding="utf-8", errors="replace").read(20000)
            except OSError:
                continue
            m = re.search(r'pseudo_type\s*=\s*"(\w+)"', head)
            if (m and m.group(1).upper() in ("US", "USPP", "PAW")) or "Ultrasoft" in head:
                uspp = True
        ew, er = self.ecut
        if ew and uspp and (er is None or er / ew <= 4 + 1e-6):
            self.add("WARN", "ultrasoft/PAW pseudopotentials with ecutrho = 4*ecutwfc; "
                             "usually 8-12*ecutwfc (see SSSP)")


def main():
    ap = argparse.ArgumentParser(description="Static checks of pw.x input files.")
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--pseudo-dir", help="override pseudo_dir of the input")
    ap.add_argument("--strict", action="store_true", help="exit 1 also on WARN")
    args = ap.parse_args()
    status = 0
    for f in args.inputs:
        if not os.path.isfile(f):
            print("%s: not found" % f)
            status = 1
            continue
        msgs = Checker(f, args.pseudo_dir).run()
        print("== %s" % f)
        if not msgs:
            print("   OK: no problems found")
        for lvl, msg, line in msgs:
            print("   %-5s %s%s" % (lvl, ("line %d: " % line) if line else "", msg))
        if any(m[0] == "ERROR" for m in msgs) or (args.strict and any(m[0] == "WARN" for m in msgs)):
            status = 1
    sys.exit(status)


if __name__ == "__main__":
    main()
