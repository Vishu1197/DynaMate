#!/usr/bin/env python3
"""
build_topology.py - patch a pdb2gmx topol.top to include cofactor and ligand.

Editing topol.top by hand is where protein-cofactor setups go wrong, because
GROMACS cares about the ORDER of includes and gives you very little help when
you get it wrong:

  * every [ atomtypes ] must appear before the first [ moleculetype ] that
    uses it - so cofactor atomtypes go directly after the forcefield include;
  * molecule topologies must come after the protein's moleculetype but BEFORE
    the water and ion includes;
  * the [ molecules ] section at the bottom must list every molecule in the
    same order as the coordinate file, or grompp assigns the wrong parameters
    to the wrong atoms without necessarily complaining.

This script does all three, and is idempotent - running it twice does not
double up the includes.

    python3 build_topology.py -p topol.top \\
        --atomtypes COF_atomtypes.itp --atomtypes ION_atomtypes.itp \\
        --atomtypes LIG_atomtypes.itp \\
        --molecule  COF_molecule.itp  --molecule  LIG_molecule.itp \\
        --molecules "Protein_chain_A:1" --molecules "COF:1" \\
        --molecules "ION:1" --molecules "LIG:1"

A backup of the original is written before anything is changed.

Part of DynaMate: https://github.com/Vishu1197/DynaMate
License: MIT
"""

import argparse
import re
import shutil
import sys
from pathlib import Path

FF_INCLUDE = re.compile(r'^\s*#include\s+".*forcefield\.itp"')
WATER_INCLUDE = re.compile(r'^\s*#include\s+".*(tip\dp|spc\w*|tip4p\w*)\.itp"',
                           re.IGNORECASE)
ION_INCLUDE = re.compile(r'^\s*#include\s+".*ions\.itp"', re.IGNORECASE)
SECTION = re.compile(r"^\s*\[\s*([A-Za-z_]+)\s*\]")

ATOMTYPE_BANNER = "; --- DynaMate: cofactor/ligand atom types ---"
MOLECULE_BANNER = "; --- DynaMate: cofactor/ligand molecule topologies ---"


def main():
    ap = argparse.ArgumentParser(
        description="Insert cofactor/ligand includes into a pdb2gmx topol.top.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("-p", "--topology", default="topol.top",
                    help="the topology to patch (default: topol.top)")
    ap.add_argument("--atomtypes", action="append", default=[], metavar="ITP",
                    help="an *_atomtypes.itp to include after the force field "
                         "(repeatable)")
    ap.add_argument("--molecule", action="append", default=[], metavar="ITP",
                    help="a *_molecule.itp to include before the water "
                         "topology (repeatable)")
    ap.add_argument("--molecules", action="append", default=[], metavar="NAME:N",
                    help="an entry for the [ molecules ] section, e.g. COF:1 "
                         "(repeatable, IN COORDINATE-FILE ORDER)")
    ap.add_argument("--no-backup", action="store_true")
    ap.add_argument("--dry-run", action="store_true",
                    help="show what would change and write nothing")
    args = ap.parse_args()

    top = Path(args.topology)
    if not top.is_file():
        sys.exit(f"Topology not found: {top}")

    for itp in args.atomtypes + args.molecule:
        if not (top.parent / itp).is_file() and not Path(itp).is_file():
            print(f"  ! {itp} does not exist yet (relative to {top.parent}). "
                  f"Including it anyway.")

    lines = top.read_text(encoding="utf-8", errors="replace").splitlines()
    original = list(lines)

    # ---------------------------------------------------------------- 1/3 --
    # Atom types: immediately after the force field include.
    to_add = [i for i in args.atomtypes
              if not any(f'#include "{i}"' in l for l in lines)]
    if to_add:
        idx = next((n for n, l in enumerate(lines) if FF_INCLUDE.match(l)), None)
        if idx is None:
            sys.exit("Could not find the forcefield.itp include in "
                     f"{top}. Is this really a pdb2gmx topology?")
        block = [""] + ([ATOMTYPE_BANNER] if ATOMTYPE_BANNER not in lines else [])
        block += [f'#include "{i}"' for i in to_add]
        lines[idx + 1:idx + 1] = block

    # ---------------------------------------------------------------- 2/3 --
    # Molecule topologies: before the water include (falling back to the ion
    # include, then to the [ system ] section).
    to_add = [i for i in args.molecule
              if not any(f'#include "{i}"' in l for l in lines)]
    if to_add:
        idx = next((n for n, l in enumerate(lines) if WATER_INCLUDE.match(l)),
                   None)
        anchor = "water topology"
        if idx is None:
            idx = next((n for n, l in enumerate(lines) if ION_INCLUDE.match(l)),
                       None)
            anchor = "ion topology"
        if idx is None:
            idx = next((n for n, l in enumerate(lines)
                        if SECTION.match(l)
                        and SECTION.match(l).group(1).lower() == "system"), None)
            anchor = "[ system ] section"
        if idx is None:
            idx = len(lines)
            anchor = "end of file"
        # Step back over any comment lines directly above the anchor so the
        # insert lands before "; Include water topology", not after it.
        while idx > 0 and lines[idx - 1].strip().startswith(";"):
            idx -= 1
        block = ([MOLECULE_BANNER] if MOLECULE_BANNER not in lines else [])
        block += [f'#include "{i}"' for i in to_add] + [""]
        lines[idx:idx] = block
        print(f"  molecule includes inserted before the {anchor}")

    # ---------------------------------------------------------------- 3/3 --
    # [ molecules ] section: rewrite it entirely if entries were supplied.
    if args.molecules:
        entries = []
        for spec in args.molecules:
            if ":" in spec:
                name, _, count = spec.rpartition(":")
            else:
                name, count = spec, "1"
            try:
                int(count)
            except ValueError:
                sys.exit(f"Bad --molecules value {spec!r} - expected NAME:COUNT")
            entries.append((name.strip(), count.strip()))

        idx = next((n for n, l in enumerate(lines)
                    if SECTION.match(l)
                    and SECTION.match(l).group(1).lower() == "molecules"), None)
        if idx is None:
            sys.exit(f"No [ molecules ] section found in {top}.")

        # Preserve any existing entries the user did not mention - typically
        # SOL and the counter-ions added later by gmx solvate / gmx genion.
        existing = []
        for l in lines[idx + 1:]:
            if SECTION.match(l):
                break
            s = l.strip()
            if not s or s.startswith(";"):
                continue
            parts = s.split()
            if len(parts) >= 2:
                existing.append((parts[0], parts[1]))
        named = {n for n, _ in entries}
        carried = [e for e in existing if e[0] not in named]

        end = idx + 1
        while end < len(lines) and not SECTION.match(lines[end]):
            end += 1

        block = ["[ molecules ]", "; Compound        #mols"]
        block += [f"{n:<18} {c}" for n, c in entries + carried]
        lines[idx:end] = block + [""]

        if carried:
            print(f"  kept existing entries: "
                  f"{', '.join(n for n, _ in carried)}")

    # ---------------------------------------------------------------- out --
    if lines == original:
        print(f"\n  {top} already has everything - nothing to do.")
        return

    if args.dry_run:
        print("\n-- dry run: the patched topology would be --\n")
        import difflib
        for d in difflib.unified_diff(original, lines,
                                      fromfile=str(top), tofile="patched",
                                      lineterm="", n=2):
            print("  " + d)
        return

    if not args.no_backup:
        bak = top.with_suffix(top.suffix + ".bak")
        n = 1
        while bak.exists():
            bak = top.with_suffix(f"{top.suffix}.bak{n}")
            n += 1
        shutil.copy2(top, bak)
        print(f"  backup written: {bak.name}")

    top.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  patched: {top}")

    print("\n" + "=" * 72)
    print("  Check this before running grompp:")
    print("=" * 72)
    for n, l in enumerate(lines, 1):
        if l.strip().startswith("#include") or SECTION.match(l):
            print(f"  {n:>5}  {l.strip()}")
    print("=" * 72)
    print("  The [ molecules ] order must match the atom order in your .gro.")


if __name__ == "__main__":
    main()
