#!/usr/bin/env python3
"""
split_receptor.py - split a receptor PDB into protein, cofactors and ions.

The usual approach is a row of greps:

    grep -E '^(ATOM)'            REC.pdb > protein.pdb
    grep -E '^(HETATM).* GDP '   REC.pdb > GDP.pdb
    grep -E '^(HETATM).* MG '    REC.pdb > MG.pdb

Those work until a residue name collides with something else on the line - an
atom named MG in a different residue, a chain identifier, a B-factor that
happens to read " MG " - at which point you silently get the wrong atoms.
This reads the PDB's fixed columns instead, and it tells you about anything it
did not extract, which is the failure that actually costs people a day: a
crystallographic water or a second ligand that vanishes without a word.

    python3 split_receptor.py -i REC.pdb --extract GDP --extract MG
        -> protein.pdb, GDP.pdb, MG.pdb

    python3 split_receptor.py -i REC.pdb --list      # just show me what's in it

Part of DynaMate: https://github.com/Vishu1197/DynaMate
License: MIT
"""

import argparse
import sys
from collections import OrderedDict
from pathlib import Path

# Residue names that are structural water, not something you want to keep.
WATERS = {"HOH", "WAT", "DOD", "TIP", "TIP3", "SOL", "H2O"}


def read_pdb(path):
    """Return (all_lines, atom_records) where records are (line, kind, resname)."""
    records = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            kind = "ATOM" if line.startswith("ATOM") else "HETATM"
            resname = line[17:20].strip().upper()
            records.append((line, kind, resname))
    return lines, records


def inventory(records):
    """{(kind, resname): [count_atoms, count_residues]} preserving order."""
    inv = OrderedDict()
    seen_res = set()
    for line, kind, resname in records:
        key = (kind, resname)
        inv.setdefault(key, [0, 0])
        inv[key][0] += 1
        # chain + residue sequence number + insertion code
        rid = (resname, line[21:22], line[22:27])
        if rid not in seen_res:
            seen_res.add(rid)
            inv[key][1] += 1
    return inv


def write_pdb(path, lines, header):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"REMARK   1 {header}\n")
        for line in lines:
            fh.write(line.rstrip() + "\n")
        fh.write("END\n")


def main():
    ap = argparse.ArgumentParser(
        description="Split a receptor PDB into protein and named "
                    "cofactor/ion components.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("-i", "--input", required=True, help="receptor PDB")
    ap.add_argument("-d", "--outdir", default=None,
                    help="output directory (default: beside the input)")
    ap.add_argument("--extract", action="append", default=[], metavar="RESNAME",
                    help="residue name to pull into its own PDB "
                         "(repeatable, e.g. --extract GDP --extract MG)")
    ap.add_argument("--protein-out", default="protein.pdb",
                    help="filename for the protein ATOM records")
    ap.add_argument("--keep-waters", action="store_true",
                    help="do not warn about crystallographic waters")
    ap.add_argument("--list", action="store_true",
                    help="show what the file contains and exit")
    args = ap.parse_args()

    src = Path(args.input)
    if not src.is_file():
        sys.exit(f"Input not found: {src}")
    outdir = Path(args.outdir) if args.outdir else src.parent
    outdir.mkdir(parents=True, exist_ok=True)

    _, records = read_pdb(src)
    if not records:
        sys.exit(f"{src}: no ATOM or HETATM records found.")

    inv = inventory(records)

    print("=" * 72)
    print(f"  {src}")
    print("=" * 72)
    print(f"  {'record':<8} {'residue':<10} {'atoms':>8} {'copies':>8}")
    print("  " + "-" * 38)
    for (kind, resname), (natoms, nres) in inv.items():
        if kind == "ATOM":
            continue
        print(f"  {kind:<8} {resname:<10} {natoms:>8} {nres:>8}")
    n_protein = sum(v[0] for k, v in inv.items() if k[0] == "ATOM")
    n_prot_res = sum(v[1] for k, v in inv.items() if k[0] == "ATOM")
    print("  " + "-" * 38)
    print(f"  {'ATOM':<8} {'(protein)':<10} {n_protein:>8} {n_prot_res:>8}")
    print("=" * 72)

    if args.list:
        return

    if not args.extract:
        print("\n  No --extract given, so only the protein will be written.")
        print("  Pick the residue names you need from the table above, e.g.:")
        het = [k[1] for k in inv if k[0] == "HETATM" and k[1] not in WATERS]
        if het:
            print("      " + " ".join(f"--extract {h}" for h in het[:4]))

    wanted = {r.strip().upper() for r in args.extract}
    missing = wanted - {k[1] for k in inv}
    if missing:
        sys.exit(f"\nRequested residue(s) not present in {src.name}: "
                 f"{', '.join(sorted(missing))}\n"
                 f"Run with --list to see what is actually there.")

    # --- write -------------------------------------------------------------
    protein_lines = [l for l, kind, rn in records
                     if kind == "ATOM" and rn not in wanted]
    written = []

    if protein_lines:
        p = outdir / args.protein_out
        write_pdb(p, protein_lines, f"protein extracted from {src.name}")
        written.append((p, len(protein_lines)))
    else:
        print("\n  ! no protein ATOM records found - is this a ligand-only file?")

    for resname in sorted(wanted):
        lines = [l for l, _, rn in records if rn == resname]
        p = outdir / f"{resname}.pdb"
        write_pdb(p, lines, f"{resname} extracted from {src.name}")
        written.append((p, len(lines)))

    print("\n  written:")
    for p, n in written:
        print(f"    {p.name:<24} {n:>7} atoms")

    # --- what did we leave behind? -----------------------------------------
    leftover = [(k[1], v[0], v[1]) for k, v in inv.items()
                if k[0] == "HETATM" and k[1] not in wanted]
    waters = [x for x in leftover if x[0] in WATERS]
    others = [x for x in leftover if x[0] not in WATERS]

    if others:
        print("\n  ! NOT EXTRACTED - these HETATM residues are in the input "
              "but not in any output:")
        for rn, natoms, nres in others:
            print(f"      {rn:<10} {natoms:>6} atoms in {nres} copies")
        print("    If any of them belong in your simulation, re-run with "
              "--extract for each.")
        print("    Otherwise they are being discarded on purpose - "
              "just be sure that is what you want.")
    if waters and not args.keep_waters:
        total = sum(x[1] for x in waters)
        print(f"\n  note: {total} crystallographic water atoms were dropped. "
              f"That is normal -")
        print("        gmx solvate will add explicit solvent later. Keep them "
              "only if a")
        print("        specific water mediates the binding you care about.")


if __name__ == "__main__":
    main()
