#!/usr/bin/env python3
"""
merge_gro.py - concatenate GROMACS .gro coordinate files into one system.

This replaces the ad-hoc `coordinate.py` / `combine_ligand.py` scripts that
most GROMACS protein-cofactor workflows end up carrying around, and it does
the two things those scripts usually get wrong: it renumbers atoms correctly
across the 5-digit wrap, and it tells you the exact [ molecules ] order the
result requires.

    python3 merge_gro.py -o complex.gro protein.gro COF.gro ION.gro
    python3 merge_gro.py -o complex_LIG.gro --box-from boxed.gro boxed.gro LIG_boxed.gro
    python3 merge_gro.py -o complex.gro --title "Protein-cofactor complex" a.gro b.gro

THE RULE THAT BREAKS EVERYONE: the order of molecules in the coordinate file
must match the order of the [ molecules ] section in topol.top, exactly. Get
it wrong and grompp either fails with an atom-count mismatch or - much worse -
succeeds while silently assigning the wrong parameters to the wrong atoms.
This script prints the required order so you can check it.

Part of DynaMate: https://github.com/Vishu1197/DynaMate
License: MIT
"""

import argparse
import sys


class GroFile:
    """A parsed GROMACS .gro file."""

    def __init__(self, path):
        self.path = path
        with open(path, encoding="utf-8", errors="replace") as fh:
            lines = fh.read().splitlines()

        if len(lines) < 3:
            sys.exit(f"{path}: too short to be a .gro file")

        self.title = lines[0].strip()
        try:
            self.natoms = int(lines[1].strip())
        except ValueError:
            sys.exit(f"{path}: line 2 should be the atom count, got {lines[1]!r}")

        body = lines[2:]
        if len(body) < self.natoms + 1:
            sys.exit(f"{path}: header says {self.natoms} atoms but only "
                     f"{len(body) - 1} coordinate lines are present")

        self.atoms = body[:self.natoms]
        self.box = body[self.natoms].rstrip()

        # Velocities are present if a coordinate line is long enough to hold
        # three more 8-char floats after the 44-char coordinate block.
        self.has_velocities = any(len(a) >= 68 for a in self.atoms)

    def residues(self):
        """Ordered list of (resnum, resname) as they appear."""
        seen = []
        last = None
        for a in self.atoms:
            key = (a[0:5].strip(), a[5:10].strip())
            if key != last:
                seen.append(key)
                last = key
        return seen

    def resnames(self):
        """Unique residue names, in order of first appearance."""
        out = []
        for _, name in self.residues():
            if name not in out:
                out.append(name)
        return out


def parse_box(box_line):
    try:
        return [float(v) for v in box_line.split()]
    except ValueError:
        sys.exit(f"Could not parse box vector line: {box_line!r}")


def main():
    ap = argparse.ArgumentParser(
        description="Merge GROMACS .gro files into a single coordinate file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("inputs", nargs="+",
                    help=".gro files to merge, IN THE ORDER they must appear "
                         "in the [ molecules ] section of topol.top")
    ap.add_argument("-o", "--output", required=True, help="output .gro file")
    ap.add_argument("--title", default=None,
                    help="title line for the merged file "
                         "(default: derived from the inputs)")
    ap.add_argument("--box-from", default=None,
                    help="take the box vectors from this file "
                         "(default: the largest box among the inputs)")
    ap.add_argument("--renumber-residues", action="store_true",
                    help="give every input a fresh residue-number range so no "
                         "two molecules share a residue number")
    args = ap.parse_args()

    if len(args.inputs) < 2:
        sys.exit("Need at least two .gro files to merge.")

    parts = [GroFile(p) for p in args.inputs]

    # --- box ---------------------------------------------------------------
    if args.box_from:
        match = [p for p in parts if p.path == args.box_from]
        if match:
            box = match[0].box
        else:
            box = GroFile(args.box_from).box
        box_src = args.box_from
    else:
        # Largest volume wins. In practice one input is the boxed receptor and
        # the rest are single molecules carrying a meaningless unit box.
        best, best_vol = None, -1.0
        for p in parts:
            v = parse_box(p.box)
            vol = v[0] * v[1] * v[2] if len(v) >= 3 else 0.0
            if vol > best_vol:
                best, best_vol = p, vol
        box, box_src = best.box, best.path

    if any(p.has_velocities for p in parts) and not all(p.has_velocities
                                                        for p in parts):
        print("  ! some inputs have velocities and some do not - velocities "
              "will be dropped", file=sys.stderr)
    keep_vel = all(p.has_velocities for p in parts)

    # --- merge -------------------------------------------------------------
    out_atoms = []
    atom_no = 0
    res_offset = 0
    summary = []

    for p in parts:
        first_atom = atom_no + 1
        local_res_max = 0
        for line in p.atoms:
            atom_no += 1
            resnum_raw = line[0:5]
            resname = line[5:10]
            atomname = line[10:15]
            coords = line[20:44]
            vel = line[44:68] if keep_vel and len(line) >= 68 else ""

            if args.renumber_residues:
                try:
                    rn = int(resnum_raw)
                except ValueError:
                    rn = 1
                local_res_max = max(local_res_max, rn)
                resnum = f"{(rn + res_offset) % 100000:5d}"
            else:
                resnum = resnum_raw

            # Atom numbers wrap at 100000 in the fixed-width .gro format.
            out_atoms.append(
                f"{resnum}{resname}{atomname}{atom_no % 100000:5d}{coords}{vel}")

        if args.renumber_residues:
            res_offset += local_res_max

        summary.append((p.path, p.natoms, first_atom, atom_no, p.resnames()))

    title = args.title or ("Merged system: " +
                           " + ".join("/".join(p.resnames()[:3])
                                      for p in parts))

    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(title[:100] + "\n")
        fh.write(f"{len(out_atoms):5d}\n")
        fh.write("\n".join(out_atoms) + "\n")
        fh.write(box + "\n")

    # --- report ------------------------------------------------------------
    print("=" * 72)
    print(f"  merged {len(parts)} files -> {args.output}")
    print(f"  box taken from : {box_src}")
    print(f"  box vectors    : {box.strip()}")
    print(f"  velocities     : {'kept' if keep_vel else 'not written'}")
    print("=" * 72)
    for path, n, a, b, names in summary:
        print(f"  {path:<28} {n:>7} atoms   #{a}-{b}   "
              f"[{', '.join(names[:4])}{' ...' if len(names) > 4 else ''}]")
    print("-" * 72)
    print(f"  {'TOTAL':<28} {len(out_atoms):>7} atoms")
    print("=" * 72)

    print("\n  The [ molecules ] section of topol.top MUST list these in this")
    print("  order. Atom order in the coordinate file and molecule order in")
    print("  the topology are matched positionally - a mismatch is silent.\n")
    print("  [ molecules ]")
    print("  ; Compound        #mols")
    for path, _, _, _, names in summary:
        if len(names) > 1:
            # Many residue names means a macromolecule: one moleculetype,
            # whose name comes from the topology, not from the residues.
            print(f"  {'<moleculetype>':<18} 1        "
                  f"; {path}: {len(names)} residue types, so this is ONE")
            print(f"  {'':<18}          ; molecule - use the name from its "
                  f"[ moleculetype ]")
            print(f"  {'':<18}          ; in topol.top, e.g. Protein_chain_A")
        else:
            print(f"  {names[0]:<18} 1        ; {path}")
    print("\n  Counts assume one copy of each. SOL and counter-ions are")
    print("  appended later by gmx solvate and gmx genion - leave them alone.")


if __name__ == "__main__":
    main()
