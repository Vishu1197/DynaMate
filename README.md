# DynaMate

**An automated GROMACS molecular dynamics workflow for protein + cofactor + metal ion + ligand systems.**

Setting up an MD simulation of a protein that carries a cofactor and a metal
ion is not hard in the sense of being conceptually difficult. It is hard in
the sense that there are about thirty places where a single wrong character
produces a topology that is *syntactically valid, silently wrong, and runs
for a hundred nanoseconds anyway*.

DynaMate is the set of scripts that grew out of doing it by hand, hitting
every one of those places, and writing down what the error actually meant.

```bash
mkdir -p ~/md/my_system && cd ~/md/my_system
cp /wherever/REC.pdb /wherever/LIG.pdb .
/path/to/DynaMate/run_all.sh
```

That is the whole interface. Two structure files in your working directory,
one command.

---

## Everything is generic

DynaMate never hardcodes a molecule name. Whatever your system actually
contains, it is renamed internally to:

| Generic name | What it is | Example |
|---|---|---|
| `Protein` | the polypeptide chain(s) | any protein |
| `COF` | the cofactor | GDP, GTP, ATP, ADP, NAD, NADP, FAD, FMN, SAM, CoA, HEM, PLP, TPP … |
| `MI` | the metal ion | MG, ZN, MN, CA, FE, CU, NI, K, NA … |
| `LIG` | the ligand under study | your docked compound |

So `COF.pdb`, `COF_molecule.itp`, `posre_COF.itp`, `rmsd_COF.xvg` mean the
same thing in every project you ever run, and the commands in the manual do
not change when the chemistry does.

You do not have to tell DynaMate what your cofactor and metal are. It reads
`REC.pdb`, works it out, and says so:

```
  split_receptor.py   REC.pdb  ->  .
==========================================================================
  Protein.pdb     5337 atoms     333 residues
  COF.pdb           40 atoms   from GDP (auto-detected), renamed to COF
                       net charge -3  (from table)
  MI.pdb             1 atoms   from MG (auto-detected), renamed to MI
                       net charge +2

  left out of the simulation:
    HOH        486 atoms   water
    GOL         12 atoms   crystallisation additive
```

Override it in `config.sh` if it guesses wrong. A system with no cofactor, no
metal, or neither is handled too — those stages just do not run.

---

## Requirements

| | |
|---|---|
| **GROMACS** ≥ 2021 | `source /usr/local/gromacs/bin/GMXRC` (2021+ for the C-rescale barostat) |
| **AmberTools** | `conda create -n ambertools -c conda-forge ambertools` |
| **Python** ≥ 3.8 | standard library only, no pip packages |
| *optional* gmx_MMPBSA | `pip install gmx_MMPBSA`, for binding free energies |

---

## Preparing the two input files

You need `REC.pdb` and `LIG.pdb`, split out of your docked complex. In UCSF
Chimera:

**REC.pdb** — protein **with** its cofactor and metal ion, **without** the
docked ligand. Select the ligand, `Actions > Atoms/Bonds > Delete`, then
`File > Save PDB`.

**LIG.pdb** — the ligand alone, **with hydrogens**:

1. `Select > Residue >` your ligand, then `Actions > Atoms/Bonds > Show only`
2. `Tools > Structure Editing > AddH` — GAFF2 needs a fully protonated
   molecule; a ligand without hydrogens gives wrong atom types and wrong
   charges, and nothing warns you
3. `Tools > Structure Editing > Minimize Structure` — a few hundred steps, to
   clean up the geometry antechamber's semi-empirical step will otherwise
   choke on
4. `File > Save PDB > Save selected atoms only`

Both files go in your working directory. Nothing else is needed.

---

## What each stage does

```
run_all.sh                        the whole thing
run_all.sh --setup                stages 1-6, stopping before production
run_all.sh --from 3               resume after fixing something
run_all.sh --list                 what is done so far
```

| Stage | Script | What happens |
|---|---|---|
| 1 | `01_prepare.sh` | split `REC.pdb` into `Protein.pdb` / `COF.pdb` / `MI.pdb`; `pdb2gmx` builds the protein topology |
| 2 | `02_parameterize.sh` | antechamber → parmchk2 → tleap → GROMACS, for COF and LIG; the metal comes from the force field |
| 3 | `03_assemble.sh` | merge coordinates, box, position restraints, patch `topol.top`, and **prove the unsolvated system is valid** |
| 4 | `04_solvate.sh` | water, then counter-ions |
| 5 | `05_minimize.sh` | steepest descents |
| 6 | `06_equilibrate.sh` | NVT 500 ps, then NPT 1 ns, both restrained |
| 7 | `07_production.sh` | unrestrained production MD, restartable from `md.cpt` |
| 8 | `08_analyze.sh` | centred trajectory, RMSD, RMSF, Rg, SASA, H-bonds, contacts, frames, windows, MM-PBSA input |

Every stage is a standalone script. If stage 3 fails, fix the problem and run
`scripts/03_assemble.sh` again — or `run_all.sh --from 3`. A finished stage
leaves a marker file and is skipped on a re-run; `FORCE=1` redoes it.

---

## The things that go wrong, and what DynaMate does about them

These are not hypothetical. Each one cost hours the first time.

**tleap writes a 0-byte topology and exits 0.**
If you forget `loadamberparams COF.frcmod`, tleap reports sixty-five
`No torsion terms for atom types: os-c3-c5-c5` errors, writes an empty
`.prmtop`, and returns success. The next command then fails with
`parmed.exceptions.FormatNotFound: Could not identify file format`, which
tells you nothing about the real cause.
→ DynaMate writes the tleap input file itself, with the `loadamberparams`
line in it, and checks the `.prmtop` is non-empty before going on.

**Fourteen duplicate-atomtype warnings, then `Fatal error: Too many warnings`.**
The cofactor and the ligand are both organic, so both GAFF tables define
`c3`, `ca`, `os`, `hn`, `o`, `n`. Include both and grompp objects to every
repeat. `-maxwarn 14` makes it run — and also silences the two warnings you
needed to read.
→ `merge_atomtypes.py` folds every table into one `gaff_atomtypes.itp`, drops
anything the force field already defines, and warns loudly if two tables
disagree about the same type. grompp then runs with zero warnings.

**`sed -n '21,447p' COF.top > COF.itp`.**
Works exactly once. Change the cofactor and those line numbers silently
produce a truncated topology — which is usually still valid, so grompp
accepts it.
→ `extract_itp.py` finds the sections by name.

**`LIG    1SOL    24741`.**
`gmx solvate` appends its `[ molecules ]` entry with a bare write. If
`topol.top` does not end in a newline, the entry lands on the previous line
and grompp reports a coordinate/topology mismatch with no hint about why.
→ every stage ensures the trailing newline, and stage 4 repairs the damage if
it happens anyway.

**`Unknown error, perhaps your text file uses wrong line endings?`**
This is what grompp says when a `posre_COF.itp` contains an `#include` of
itself — the result of pasting the `#ifdef POSRES` block into the wrong file.
→ the block is written into `COF_molecule.itp` automatically, where it
belongs.

**`gmx editconf -translate 8.156 8.255 5.515`.**
Boxing the receptor first and then moving the ligand to match means copying
three numbers out of editconf's output by hand. They are right for that one
run. Re-run with a different box and the ligand ends up in the solvent —
valid file, no error, flat RMSD, wasted week.
→ DynaMate merges the ligand **before** boxing, in the frame it was docked
in, and `combine.py` refuses to continue if the ligand is not actually
touching the protein.

**`gmx genion` group 18.**
Group numbers depend on what is in the system and change the moment anything
does.
→ every interactive selection is made by **name**, not number.

**Production MD with `define = -DPOSRES` still set.**
The run completes normally, the files look right, and the protein never
moved.
→ stage 7 refuses to start if `md.mdp` still has it.

**`-dump 25` is 25 picoseconds, not 25 nanoseconds.**
A frame at 25 ps certainly exists, so nothing complains.
→ stage 8 converts nanoseconds to picoseconds and prints both.

**A 1.33 nm bond through the middle of the protein.**
pdb2gmx bonds consecutive residues whatever the distance, so a disordered
loop with no density becomes a rubber band. It says `Long Bond (1151-1153 =
1.3316 nm)` once and moves on; grompp mentions it much later as
`largest distance between excluded atoms is 1.459 nm`.
→ stage 1 catches it and explains what it means.

**`System has non-zero total charge: 2.997002`.**
The `0.003` is not an error and no number of ions will remove it: AM1-BCC
partial charges are rounded per atom, so a −3 cofactor sums to −2.999.
→ reported once, explained, and not treated as a problem.

---

## Configuration

`config.sh` is the only file you normally edit, and most of it has a sensible
default. The ones worth reading before your first run:

```bash
COF_RESNAME="auto"        # or GDP, ATP, NAD ... or "none"
COF_CHARGE="auto"         # looked up for ~40 common cofactors
MI_RESNAME="auto"         # or MG, ZN, MN ... or "none"
LIG_CHARGE=0              # COUNT THIS YOURSELF

FF_NAME="amber99sb-ildn"  # by name - menu numbers move between versions
WATER_MODEL="tip3p"

MI_SOURCE="forcefield"    # take the ion from the force field's ions.itp
                          # "amber" for the Li/Merz set via tleap

MDRUN_EXTRA=""            # "-nb gpu -pme gpu -ntmpi 1" for one GPU
```

`LIG_CHARGE` is the one number nothing can check for you. Neutral drug-like
molecule 0; one carboxylate −1; one protonated amine +1. Get it wrong and
every partial charge on the ligand is wrong, in a topology that looks fine.

Simulation length is in `mdp/md.mdp`:

```
nsteps = ns × 500000        (with dt = 0.002 ps)
 50 ns  ->  25000000
100 ns  ->  50000000
500 ns  -> 250000000
```

---

## Output

```
my_system/
├── REC.pdb  LIG.pdb              what you provided
├── Protein.pdb  COF.pdb  MI.pdb  the split receptor
├── topol.top                     the assembled topology
├── gaff_atomtypes.itp            one de-duplicated atom type table
├── COF_molecule.itp  MI_molecule.itp  LIG_molecule.itp
├── posre_COF.itp  posre_LIG.itp
├── complex.gro  complex_LIG.gro  boxed.gro  solvated_ions.gro
├── em.gro  nvt.gro  npt.gro  md.gro
├── md.xtc  md.edr  md.log  md.cpt
├── dynamate.log  dynamate.env    what ran, and what was detected
└── analysis/
    ├── md_center.xtc  index.ndx  start.pdb
    ├── rmsd_backbone.xvg  rmsd_COF.xvg  rmsd_LIG.xvg
    ├── rmsf_residue.xvg  gyrate.xvg  sasa.xvg
    ├── hbond_protein_LIG.xvg  mindist_protein_LIG.xvg
    ├── mindist_MI_COF.xvg     is the metal still coordinated?
    ├── frame_0ns.pdb  frame_25ns.pdb  …
    ├── traj_0_10ns.xtc  traj_90_100ns.xtc  …
    └── mmpbsa.in  run_mmpbsa.sh
```

---

## Documentation

- **`docs/WORKFLOW.md`** — the whole procedure as explicit commands, for
  running by hand or for understanding what the scripts do
- **`docs/TROUBLESHOOTING.md`** — every error message this workflow can
  produce, what it actually means, and how to fix it

---

## A note on the physics

DynaMate automates the mechanics. It does not make the scientific decisions
for you, and three of them are yours alone:

- **The cofactor's net charge.** Auto-detection covers common cofactors, but
  the protonation state in *your* structure at *your* pH is your call.
- **The ligand's protonation state.** Both the charge and which atoms carry
  the hydrogens.
- **The metal model.** DynaMate uses a non-bonded point charge, which is the
  standard treatment and what nearly all published Mg²⁺-nucleotide work does.
  It does not enforce coordination geometry, and the ion can leave its site
  during a long run. That is a real result, not a bug — but if you need
  enforced coordination, you need a bonded or 12-6-4 model, and you need to
  build it deliberately.

Check `analysis/mindist_MI_COF.xvg` before you believe anything about the
metal site.

---

## Citation

If DynaMate is useful in your work, please cite the repository and the tools
it orchestrates — GROMACS, AmberTools/GAFF2, and gmx_MMPBSA if you use it.

```
Chanda, V. DynaMate: an automated GROMACS molecular dynamics workflow for
protein-cofactor-metal-ligand systems. v2.0 (2026).
https://github.com/Vishu1197/DynaMate
```

**Vishal Chanda** — Doctoral researcher, REVA University, Bengaluru
ORCID [0000-0002-1646-6329](https://orcid.org/0000-0002-1646-6329) ·
[Google Scholar](https://scholar.google.com/citations?user=gdA8TAIAAAAJ) ·
[ResearchGate](https://www.researchgate.net/profile/Vishal-Chanda)

## License

MIT — see `LICENSE`.
