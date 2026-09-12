#!/usr/bin/env python3
"""
extract_itp.py - split an Amber-converted GROMACS .top into includable .itp files.

When you parameterise a cofactor or ligand with antechamber/tleap and convert
the result with amb2gro_top_gro.py (or ACPYPE), you get a standalone .top that
cannot simply be #included into your protein topology: it carries its own
[ defaults ], [ system ] and [ molecules ] sections, which collide with the
ones pdb2gmx already wrote.

The usual fix is to cut it up with hardcoded line numbers:

    sed -n '21,447p' GDP.top > GDP.itp

That works exactly once. Change the ligand, change the AmberTools version,
change anything, and those numbers silently produce a truncated topology -
which grompp will often accept, because the result is still syntactically
valid. This script finds the sections by name instead.

    python3 extract_itp.py -i COF.top -n COF
        -> COF_atomtypes.itp   ([ atomtypes ] only)
        -> COF_molecule.itp    ([ moleculetype ] onwards, minus [ system ])

    python3 extract_itp.py -i ION.top -n ION --atomtypes-only

Part of DynaMate: https://github.com/Vishu1197/DynaMate
License: MIT
"""

import argparse
import re
import sys
from pathlib import Path

SECTION = re.compile(r"^\s*\[\s*([A-Za-z_]+)\s*\]")

# Sections that must NOT survive into an includable .itp, because pdb2gmx has
# already provided them and GROMACS refuses duplicates.
DROP_SECTIONS = {"defaults", "system", "molecules"}


def parse_sections(path):
    """Return [(section_name, [lines including the header])] in file order."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()

    blocks, current, name = [], [], None
    for line in lines:
        m = SECTION.match(line)
        if m:
            if name is not None or current:
                blocks.append((name, current))
            name, current = m.group(1).lower(), [line]
        else:
            current.append(line)
    if name is not None or current:
        blocks.append((name, current))
    return blocks


def write_block(path, header, blocks):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(header)
        for _, lines in blocks:
            fh.write("\n".join(lines).rstrip() + "\n\n")


def main():
    ap = argparse.ArgumentParser(
        description="Split an Amber-converted .top into atomtypes and "
                    "molecule .itp files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("-i", "--input", required=True,
                    help="the .top produced by amb2gro_top_gro.py or ACPYPE")
    ap.add_argument("-n", "--name", default=None,
                    help="output prefix (default: the input file's stem)")
    ap.add_argument("-d", "--outdir", default=None,
                    help="output directory (default: beside the input)")
    ap.add_argument("--atomtypes-only", action="store_true",
                    help="write only <name>_atomtypes.itp. Use this for a "
                         "metal ion whose molecule definition already exists "
                         "in your force field's ions.itp")
    ap.add_argument("--posres", default=None,
                    help="append an '#ifdef POSRES / #include \"<file>\"' "
                         "block to the molecule .itp")
    args = ap.parse_args()

    src = Path(args.input)
    if not src.is_file():
        sys.exit(f"Input not found: {src}")

    name = args.name or src.stem
    outdir = Path(args.outdir) if args.outdir else src.parent
    outdir.mkdir(parents=True, exist_ok=True)

    blocks = parse_sections(src)
    found = [b[0] for b in blocks if b[0]]
    if not found:
        sys.exit(f"{src}: no [ sections ] found - is this really a .top file?")

    atomtypes = [b for b in blocks if b[0] == "atomtypes"]
    if not atomtypes:
        print(f"  ! {src.name} has no [ atomtypes ] section.")
        print("    That is normal for some ion topologies - the force field "
              "already defines them.")

    dropped = sorted({b[0] for b in blocks if b[0] in DROP_SECTIONS})

    # Everything from the first [ moleculetype ] onwards, minus the sections
    # that would collide with the protein topology.
    mol_start = next((i for i, b in enumerate(blocks)
                      if b[0] == "moleculetype"), None)
    molecule = []
    if mol_start is not None:
        molecule = [b for b in blocks[mol_start:] if b[0] not in DROP_SECTIONS]

    banner = (f"; Extracted from {src.name} by extract_itp.py (DynaMate)\n"
              f"; Do not edit by hand - re-run the extraction instead.\n\n")

    written = []

    if atomtypes:
        p = outdir / f"{name}_atomtypes.itp"
        write_block(p, banner, atomtypes)
        written.append(p)

    if not args.atomtypes_only:
        if not molecule:
            sys.exit(f"{src}: no [ moleculetype ] section found - "
                     f"nothing to write as a molecule topology.")
        p = outdir / f"{name}_molecule.itp"
        write_block(p, banner, molecule)
        if args.posres:
            with open(p, "a", encoding="utf-8") as fh:
                fh.write("#ifdef POSRES\n")
                fh.write(f'#include "{args.posres}"\n')
                fh.write("#endif\n")
        written.append(p)

    # --- report ------------------------------------------------------------
    print("=" * 72)
    print(f"  source   : {src}")
    print(f"  sections : {', '.join(found)}")
    if dropped:
        print(f"  dropped  : {', '.join(dropped)}  "
              f"(pdb2gmx already provides these)")
    print("=" * 72)
    for p in written:
        secs = [b[0] for b in parse_sections(p) if b[0]]
        print(f"  {p.name:<28} [{'] ['.join(secs)}]")
    if args.posres:
        print(f"  {'':<28} + #ifdef POSRES -> {args.posres}")
    print("=" * 72)

    if atomtypes and not args.atomtypes_only:
        print("\n  Include them in topol.top like this:\n")
        print("    ; right after the forcefield include")
        print(f'    #include "{name}_atomtypes.itp"')
        print("\n    ; after the protein moleculetype, "
              "BEFORE the water topology")
        print(f'    #include "{name}_molecule.itp"')
        print("\n  Never #include the original .top - it will collide.")


if __name__ == "__main__":
    main()
