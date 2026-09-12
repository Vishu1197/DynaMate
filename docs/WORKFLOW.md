# The full protocol, step by step

This is what DynaMate's scripts do, written out as commands you could type
yourself. Read it if you want to understand the pipeline, adapt it to a system
it doesn't quite fit, or debug a stage that failed.

Throughout, `COF` is your cofactor's residue name (GDP, ATP, NAD…), `ION` is
your metal ion (MG, ZN, MN…), and `LIG` is your docked ligand.

**Contents**

1. [Separating receptor from ligand](#1-separating-receptor-from-ligand)
2. [Splitting the receptor](#2-splitting-the-receptor)
3. [Protein topology](#3-protein-topology)
4. [Parameterising the cofactor](#4-parameterising-the-cofactor)
5. [Parameterising the metal ion](#5-parameterising-the-metal-ion)
6. [Parameterising the ligand](#6-parameterising-the-ligand)
7. [Assembling coordinates](#7-assembling-coordinates)
8. [Assembling the topology](#8-assembling-the-topology)
9. [Validation before solvation](#9-validation-before-solvation)
10. [Solvation and ions](#10-solvation-and-ions)
11. [Energy minimisation](#11-energy-minimisation)
12. [Position restraints](#12-position-restraints)
13. [NVT equilibration](#13-nvt-equilibration)
14. [NPT equilibration](#14-npt-equilibration)
15. [Production MD](#15-production-md)
16. [Post-processing and analysis](#16-post-processing-and-analysis)
17. [MM-PBSA binding free energy](#17-mm-pbsa-binding-free-energy)

---

## 1. Separating receptor from ligand

**This is the one step that needs a human.** Everything after it is scripted.

You have a docked complex: a protein that crystallised with a cofactor and a
metal ion, plus a drug candidate placed by docking. GROMACS needs these as
separate files because they take completely different parameterisation routes.

In **UCSF Chimera** (or ChimeraX, or PyMOL):

1. Open the complete protein–ligand complex.
2. Select the protein, the cofactor and the metal ion together.
   Save as **`REC.pdb`**.
3. Select the docked ligand alone.
   - `Tools → Structure Editing → AddH` — antechamber needs every hydrogen.
   - `Tools → Structure Editing → Minimize Structure` — a docking pose can have
     slightly strained geometry that upsets the charge fit.
   - Save as **`LIG.pdb`**.
4. Put both files in your working directory.

> **Why minimise the ligand but not the receptor?**
> The receptor's geometry comes from experiment and is treated by a
> well-validated protein force field. The ligand is about to have its charges
> fitted quantum-mechanically, and that fit is sensitive to geometry.

**Keep or discard crystallographic waters?** Discard them by default —
`gmx solvate` will add explicit solvent shortly. Keep a water only if it
specifically mediates the interaction you care about, in which case treat it as
another component with its own residue name.

---

## 2. Splitting the receptor

First, see what you actually have:

```bash
python3 scripts/split_receptor.py -i REC.pdb --list
```

```
  record   residue       atoms   copies
  --------------------------------------
  HETATM   COF              40        1
  HETATM   ION               1        1
  HETATM   HOH             112      112
  --------------------------------------
  ATOM     (protein)      5337      333
```

Then split:

```bash
python3 scripts/split_receptor.py -i REC.pdb --extract COF --extract ION
```

Producing `protein.pdb`, `COF.pdb`, `ION.pdb`.

<details>
<summary><b>Doing it with grep instead</b></summary>

The traditional approach:

```bash
grep -E '^ATOM'              REC.pdb > protein.pdb
grep -E '^HETATM.* COF '     REC.pdb > COF.pdb
grep -E '^HETATM.* ION '     REC.pdb > ION.pdb
```

This works, but matches anywhere on the line rather than in the residue-name
columns (18–20). An atom named MG inside a different residue, or a chain ID or
B-factor that happens to read ` MG `, silently ends up in the wrong file. It
also says nothing about what it left behind — a second ligand or a metal you
forgot about disappears without a word. The script reads the columns and
reports leftovers.

</details>

---

## 3. Protein topology

```bash
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -p topol.top -i posre.itp -ignh
```

You'll be asked to pick a force field and a water model. **The menu numbers
change between GROMACS versions** — read the list rather than trusting a number
from a tutorial.

For a protein–ligand system parameterised with GAFF2, **AMBER99SB-ILDN** with
**TIP3P** is the conventional pairing and what most of the literature uses.

`-ignh` discards existing hydrogens and rebuilds them according to the force
field's own naming, which is almost always what you want for a crystal
structure.

Note what it reports: residue count, atom count, and **total charge**. You'll
need that charge later.

<details>
<summary><b>If pdb2gmx fails</b></summary>

| Message | Cause |
|---|---|
| `Residue 'XXX' not found in residue topology database` | A non-standard residue is still in `protein.pdb`. Either extract it as a cofactor, or provide parameters for it |
| `Atom X in residue Y not found` | Missing backbone atoms — an unresolved loop. Model it in first, or truncate the chain |
| `Chain identifier 'A' was used in two non-sequential blocks` | Non-contiguous chains. Renumber, or split and recombine |

</details>

---

## 4. Parameterising the cofactor

GROMACS has no parameters for GDP, ATP or NAD. GAFF2 provides them, via
AmberTools.

```bash
conda activate ambertools
```

**Charges.** AM1-BCC semi-empirical charge fitting:

```bash
antechamber -i COF.pdb -fi pdb -o COF.mol2 -fo mol2 \
            -at gaff2 -c bcc -nc -3 -rn COF -pf y
```

`-nc -3` is the **net formal charge**. This is not a detail:

| Cofactor | Charge at pH 7 |
|---|---|
| GDP, ADP | −3 |
| GTP, ATP | −4 |
| NAD⁺ | −1 |
| NADH | −2 |
| FAD | −2 |

Get it wrong and antechamber succeeds, produces a valid topology, and every
partial charge in the molecule is wrong.

**Missing parameters.** GAFF2 won't have every bonded term:

```bash
parmchk2 -i COF.mol2 -f mol2 -o COF.frcmod -s 2
```

Read `COF.frcmod`. Terms marked `ATTN, need revision` are guesses. For a
common cofactor there should be few or none; a page of them means GAFF2 doesn't
really cover your molecule and you should consider published parameters
instead.

**Amber topology:**

```bash
tleap -f leap_COF.in
```

```
source leaprc.gaff2
loadamberparams COF.frcmod
COF = loadmol2 COF.mol2
check COF
charge COF
saveamberparm COF COF.prmtop COF.inpcrd
quit
```

`check` reports missing parameters and close contacts; `charge` prints the total
charge, which should match your `-nc`.

**Convert to GROMACS:**

```bash
amb2gro_top_gro.py -p COF.prmtop -c COF.inpcrd \
                   -t COF.top -g COF.gro -b COF_converted.pdb
```

(`acpype` does the same job if you prefer it.)

**Extract includable topology.** `COF.top` is standalone — it has its own
`[ defaults ]`, `[ system ]` and `[ molecules ]` sections that collide with
what pdb2gmx wrote. Split it:

```bash
python3 scripts/extract_itp.py -i COF.top -n COF
```

Producing `COF_atomtypes.itp` and `COF_molecule.itp`.

> **Why not `sed`?**
> The traditional command is `sed -n '21,447p' COF.top > COF.itp`. Those line
> numbers are specific to one molecule, one AmberTools version, one run. Change
> anything and you get a **truncated topology that is still syntactically
> valid** — grompp accepts it, and you simulate a molecule missing half its
> dihedrals. `extract_itp.py` finds sections by name.

---

## 5. Parameterising the metal ion

A monatomic ion needs no GAFF treatment — just Amber's ion parameter set.

```bash
tleap -f leap_ION.in
```

```
source leaprc.gaff2
loadamberparams frcmod.ions234lm_126_tip3p
loadOff atomic_ions.lib
ION = loadpdb ION.pdb
check ION
charge ION
saveamberparm ION ION.prmtop ION.inpcrd
quit
```

`frcmod.ions234lm_126_tip3p` covers 2+/3+/4+ metals with the 12-6 Lennard-Jones
set for TIP3P water. **Match the water model** — using TIP3P ion parameters with
TIP4P water is a real, if subtle, error.

Convert, then take **atom types only**:

```bash
amb2gro_top_gro.py -p ION.prmtop -c ION.inpcrd -t ION.top -g ION.gro -b ION_converted.pdb
python3 scripts/extract_itp.py -i ION.top -n ION --atomtypes-only
```

> [!IMPORTANT]
> **Do not include a molecule topology for the metal.** Your force field's
> `ions.itp` already defines MG, ZN, CA and the rest. Including both gives a
> duplicate `[ moleculetype ]`, and the resulting error message is genuinely
> hard to interpret. You need the atom types (so the ion's Lennard-Jones
> parameters are the Amber ones) and nothing else.

---

## 6. Parameterising the ligand

Identical to the cofactor:

```bash
antechamber -i LIG.pdb -fi pdb -o LIG.mol2 -fo mol2 -at gaff2 -c bcc -nc 0 -rn LIG -pf y
parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod -s 2
tleap -f leap_LIG.in
amb2gro_top_gro.py -p LIG.prmtop -c LIG.inpcrd -t LIG.top -g LIG.gro -b LIG_converted.pdb
python3 scripts/extract_itp.py -i LIG.top -n LIG
```

Most neutral drug-like molecules are `-nc 0`. Carboxylates are −1, amines
protonated at pH 7 are +1. Decide deliberately.

```bash
conda deactivate    # back to your GROMACS environment
```

---

## 7. Assembling coordinates

```bash
python3 scripts/merge_gro.py -o complex.gro \
        --box-from protein_processed.gro \
        protein_processed.gro COF.gro ION.gro LIG.gro
```

> [!WARNING]
> **The order of files here is the order of `[ molecules ]` in `topol.top`.**
> GROMACS matches coordinates to topology **positionally**, not by name. Get
> the order wrong and grompp may accept it while assigning your cofactor's
> parameters to your ligand's atoms.
>
> `merge_gro.py` prints the required order after every run. Compare it against
> the bottom of `topol.top` before you go further.

<details>
<summary><b>If the ligand ends up in the wrong place</b></summary>

`tleap` sometimes recentres a molecule. If your ligand's coordinates no longer
match its docked pose, translate it back:

```bash
gmx editconf -f LIG.gro -o LIG_placed.gro -translate <dx> <dy> <dz>
```

Get the offset by comparing a reference atom's coordinates in `LIG.pdb`
(Ångström) with `LIG.gro` (nanometres — divide by 10).

Better: check `LIG_converted.pdb` against the original `LIG.pdb` in a viewer
before merging, and you'll catch it immediately.

</details>

---

## 8. Assembling the topology

Three edits to `topol.top`, in three specific places.

**Atom types — directly after the force-field include.** Every atom type must
be defined before the first `[ moleculetype ]` that uses it.

```
#include "amber99sb-ildn.ff/forcefield.itp"

; cofactor / ion / ligand atom types
#include "COF_atomtypes.itp"
#include "ION_atomtypes.itp"
#include "LIG_atomtypes.itp"
```

**Molecule topologies — after the protein, before the water.**

```
#ifdef POSRES
#include "posre.itp"
#endif

; cofactor / ligand molecule topologies
#include "COF_molecule.itp"
#include "LIG_molecule.itp"

; Include water topology
#include "amber99sb-ildn.ff/tip3p.itp"
```

Note there is no `ION_molecule.itp` — see step 5.

**The molecules section — at the bottom, in coordinate-file order.**

```
[ molecules ]
; Compound        #mols
Protein_chain_A     1
COF                 1
ION                 1
LIG                 1
```

All three at once, with a backup and idempotency:

```bash
python3 scripts/build_topology.py -p topol.top \
    --atomtypes COF_atomtypes.itp --atomtypes ION_atomtypes.itp \
    --atomtypes LIG_atomtypes.itp \
    --molecule COF_molecule.itp --molecule LIG_molecule.itp \
    --molecules "Protein_chain_A:1" --molecules "COF:1" \
    --molecules "ION:1" --molecules "LIG:1"
```

Add `--dry-run` to see the diff first.

---

## 9. Validation before solvation

Test the dry system before adding 25,000 water molecules to it:

```bash
gmx grompp -f mdp/em.mdp -c complex.gro -p topol.top -o validate.tpr -maxwarn 2
```

If this passes, your topology and coordinates agree. If it fails, fixing it now
takes seconds; fixing it after solvation means redoing everything downstream.

Expect to see `System has non-zero total charge` with a small value — the
rounding residual of the fitted partial charges. Anything under ~0.01 is fine.

---

## 10. Solvation and ions

```bash
gmx editconf -f complex.gro -o boxed.gro -c -d 1.0 -bt dodecahedron
```

A dodecahedron holds the same clearance in ~29% less volume than a cube, which
is ~29% less water to simulate. `-d 1.0` is the minimum that stops the protein
from interacting with its own periodic image at a 1.0 nm cutoff.

```bash
gmx solvate -cp boxed.gro -cs spc216.gro -p topol.top -o solvated.gro
```

Use `-cs spc216.gro` rather than an absolute path — GROMACS resolves it through
its own share directory, so the command works on any machine.

```bash
gmx grompp -f mdp/ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2
echo SOL | gmx genion -s ions.tpr -o solvated_ions.gro -p topol.top \
                      -pname NA -nname CL -neutral
```

**Answer `SOL` when asked which group to take ions from.** Any other answer
replaces protein or cofactor atoms with ions.

Add `-conc 0.15` for physiological ionic strength on top of neutralisation.

Then read the result:

```bash
grep -A10 "\[ molecules \]" topol.top
```

Every component should be there exactly once, plus SOL and your counter-ions.

---

## 11. Energy minimisation

```bash
gmx grompp -f mdp/em.mdp -c solvated_ions.gro -p topol.top -o em.tpr -maxwarn 2
gmx mdrun -v -deffnm em
```

Success looks like:

```
Steepest Descents converged to Fmax < 1000 in NNNN steps
Potential Energy  = -1.2e+06
Maximum force     =  9.5e+02
```

A **positive** potential energy or a **failure to converge** means a bad
contact. Usual suspects: the docked pose overlaps the protein, or the cofactor
came back from `tleap` in a different coordinate frame. Open `em.gro` in a
viewer.

---

## 12. Position restraints

`pdb2gmx` wrote `posre.itp` for the protein. The cofactor and ligand need their
own, or they'll drift out of the binding site during the equilibration that is
supposed to be holding everything still.

```bash
echo 0 | gmx genrestr -f COF.gro -o posre_COF.itp -fc 1000 1000 1000
echo 0 | gmx genrestr -f LIG.gro -o posre_LIG.itp -fc 1000 1000 1000
```

Append to each molecule topology:

```bash
cat >> COF_molecule.itp <<'EOF'

#ifdef POSRES
#include "posre_COF.itp"
#endif
EOF
```

…and the same for `LIG_molecule.itp`.

**Leave the metal ion unrestrained.** It's a single atom whose coordination
geometry should be free to relax with the protein around it. Restraining it
fights the very relaxation you're equilibrating for.

---

## 13. NVT equilibration

Constant Number, Volume, Temperature — bring the system to temperature while
the solute is pinned, so the solvent orders itself around it.

```bash
gmx grompp -f mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2
gmx mdrun -deffnm nvt
```

**`-r` is required** whenever `define = -DPOSRES`. It supplies the reference
coordinates the restraints pull toward. Omit it and grompp errors out.

500 ps at 300 K. Check the temperature settled:

```bash
echo Temperature | gmx energy -f nvt.edr -o nvt_temperature.xvg
```

It should reach the target within a few tens of ps and then fluctuate by a few K.

---

## 14. NPT equilibration

Constant Number, Pressure, Temperature — let the box find the right density.

```bash
gmx grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt \
           -p topol.top -o npt.tpr -maxwarn 2
gmx mdrun -deffnm npt
```

`-t nvt.cpt` carries the velocities over, so this genuinely continues from NVT
rather than restarting.

1 ns at 300 K and 1 bar. Check four things:

```bash
echo Temperature | gmx energy -f npt.edr -o npt_temperature.xvg
echo Pressure    | gmx energy -f npt.edr -o npt_pressure.xvg
echo Density     | gmx energy -f npt.edr -o npt_density.xvg
echo Volume      | gmx energy -f npt.edr -o npt_volume.xvg
```

- **Temperature** — at target, fluctuating a few K
- **Pressure** — fluctuates by hundreds of bar. This is normal and physical;
  only the average is meaningful. Don't try to fix it.
- **Density** — should be *flat* by the end, around 1000 kg/m³
- **Volume** — should stop drifting

A density still climbing at the end means the box hasn't equilibrated. Extend
NPT rather than starting production on an unequilibrated system.

---

## 15. Production MD

```bash
gmx grompp -f mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 1
gmx mdrun -s md.tpr -deffnm md -cpt 15 -cpo md.cpt
```

Add hardware flags as appropriate: `-nb gpu -pme gpu` for a GPU, `-nt 12` to
cap CPU threads.

**No `-r`**, because there are no position restraints. Check `md.mdp` has
`define =` with nothing after it. A stray `-DPOSRES` here simulates a frozen
protein for the entire run, which you will not notice until analysis.

**Length**: `nsteps × dt`. At `dt = 0.002 ps`, `nsteps = ns × 500000`.
`nsteps = 50000000` is 100 ns.

### If it stops

Power cut, reboot, crash, `Ctrl-C`:

```bash
ls -lh md.cpt
gmx mdrun -s md.tpr -deffnm md -cpi md.cpt -cpt 15 -cpo md.cpt
```

> [!CAUTION]
> **Do not regenerate `md.tpr`.** It contains the run's state definition; a new
> one restarts from zero and discards everything computed so far. `-cpi md.cpt`
> is what makes the run continue toward its original endpoint.

---

## 16. Post-processing and analysis

**Fix periodic boundaries first.** Every analysis below is wrong on a raw
trajectory, because molecules jump across box edges and RMSD reads those jumps
as enormous motion.

```bash
printf 'Protein\nSystem\n' | \
gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol -ur compact
```

Then:

```bash
# Protein backbone RMSD - is the fold stable?
printf 'Backbone\nBackbone\n' | gmx rms -s md.tpr -f md_center.xtc -o rmsd.xvg -tu ns

# Ligand RMSD, fitted on the protein - DID THE DRUG STAY BOUND?
printf 'Backbone\nLIG\n' | gmx rms -s md.tpr -f md_center.xtc -o rmsd_ligand.xvg -tu ns

# Per-residue fluctuation - which loops are flexible?
echo Protein | gmx rmsf -s md.tpr -f md_center.xtc -o rmsf.xvg -res

# Compactness - did the protein unfold?
echo Protein | gmx gyrate -s md.tpr -f md_center.xtc -o gyrate.xvg

# Hydrogen bonds
printf 'Protein\nLIG\n' | gmx hbond -s md.tpr -f md_center.xtc -num hb.xvg -tu ns

# Metal coordination - did the ion stay in its site?
printf 'Protein\nION\n' | gmx mindist -s md.tpr -f md_center.xtc -od mindist_ion.xvg
```

The **ligand RMSD fitted on the protein** is the single most informative plot
for a docking follow-up. Flat and low means the pose held. A step change means
it moved to a different sub-pocket. A steady climb means it left.

> **Zero protein–ligand hydrogen bonds isn't a bug** if your ligand is a pure
> hydrocarbon — it has no donors or acceptors. Judge binding by ligand RMSD and
> contact distances instead, and use hydrogen-bond analysis for the
> cofactor–protein interaction where it's informative.

### Extracting frames and windows

`-dump` and `-b`/`-e` take **picoseconds**, not nanoseconds:

```bash
echo System | gmx trjconv -s md.tpr -f md_center.xtc -o frame_50ns.pdb -dump 50000
echo System | gmx trjconv -s md.tpr -f md_center.xtc -o traj_90_100ns.xtc -b 90000 -e 100000
```

| Time | Value to pass |
|---|---|
| 0 ns | 0 |
| 25 ns | 25000 |
| 50 ns | 50000 |
| 100 ns | 100000 |

`-dump 50` gives you the frame at 50 **picoseconds**, which is a very easy
mistake to make and a confusing one to diagnose.

---

## 17. MM-PBSA binding free energy

```bash
pip install gmx_MMPBSA
gmx_MMPBSA --help
```

Build an index file with the two groups you want:

```bash
gmx make_ndx -f md.tpr -o index.ndx
# note the group NUMBERS for your receptor and ligand, then: q
```

`mmpbsa.in`:

```
&general
startframe=1,
endframe=100,
interval=1,
verbose=1,
/
&gb
igb=5, saltcon=0.150,
/
```

> [!IMPORTANT]
> **The exact `gmx_MMPBSA` invocation depends on the version installed**, which
> is why this is documented rather than scripted. Check `gmx_MMPBSA --help`
> against your version before running it.
>
> Two things that are wrong in most copied-around workflows:
>
> - **`-c` takes a `.tpr`, not a `.pdb`.** Passing an extracted frame PDB as
>   `-c` doesn't give it the topology it needs.
> - **`-cs`/`-ct` must match.** A structure file and a trajectory that came
>   from different systems will produce numbers that look plausible and mean
>   nothing.
>
> A typical call over a trajectory window:
>
> ```bash
> gmx_MMPBSA -O -i mmpbsa.in \
>            -cs md.tpr -ct traj_90_100ns.xtc -ci index.ndx \
>            -cg <receptor_group> <ligand_group> \
>            -cp topol.top \
>            -o FINAL_RESULTS.dat -eo FINAL_RESULTS.csv
> ```
>
> Replace `<receptor_group>` and `<ligand_group>` with the numbers from
> `make_ndx`.

**On interpreting the result:** MM-PBSA and MM-GBSA binding energies are useful
for *ranking* compounds against each other in the same system. They are not
reliable absolute binding free energies — they neglect explicit solvent and
treat entropy crudely. Report them as relative, and say which method and
parameters you used.

---

## Appendix: what each file is for

| File | Written by | Purpose |
|---|---|---|
| `protein.pdb` | 01 | Protein only, for pdb2gmx |
| `COF.pdb` / `ION.pdb` | 01 | Cofactor and metal, for AmberTools |
| `COF.mol2` | antechamber | Structure with AM1-BCC partial charges |
| `COF.frcmod` | parmchk2 | Bonded parameters GAFF2 was missing |
| `COF.prmtop` / `.inpcrd` | tleap | Amber-format topology and coordinates |
| `COF.top` / `.gro` | amb2gro | Same, converted to GROMACS — **not includable as-is** |
| `COF_atomtypes.itp` | extract_itp | `[ atomtypes ]` only, goes after the force field |
| `COF_molecule.itp` | extract_itp | `[ moleculetype ]` onwards, goes before the water |
| `topol.top` | pdb2gmx + patched | The master topology |
| `complex.gro` | merge_gro | All components, in `[ molecules ]` order |
| `solvated_ions.gro` | genion | Ready to minimise |
| `em.gro` | mdrun | Minimised |
| `nvt.gro` / `npt.gro` | mdrun | Equilibrated |
| `md.tpr` | grompp | The production run definition — **never regenerate mid-run** |
| `md.xtc` | mdrun | Compressed production trajectory |
| `md.cpt` | mdrun | Restart checkpoint |
