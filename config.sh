#!/usr/bin/env bash
# =============================================================================
#  DynaMate configuration
#  ---------------------------------------------------------------------------
#  This is the only file you need to edit. Every stage script reads it.
#  Copy your project's structures into $WORKDIR, set the names below, and run
#  the stages in order.
# =============================================================================

# --- project -----------------------------------------------------------------

# Working directory. Everything is created here. Use an absolute path.
# CHANGEME is a sentinel - the scripts refuse to run until you replace it.
WORKDIR="CHANGEME"

# A short name for the system, used in titles and log lines.
SYSTEM_NAME="protein-cofactor-ligand complex"


# --- input structures --------------------------------------------------------
# Put these in $WORKDIR before you start.
#
#   RECEPTOR_PDB : protein WITH its crystallographic cofactor and metal ion
#                  still attached. Separate it from the docked ligand in
#                  Chimera/PyMOL first - see docs/WORKFLOW.md step 1.
#   LIGAND_PDB   : the docked ligand alone, hydrogens added and minimised.

RECEPTOR_PDB="REC.pdb"
LIGAND_PDB="LIG.pdb"


# --- what is in the receptor -------------------------------------------------
# Residue names EXACTLY as they appear in columns 18-20 of your RECEPTOR_PDB.
# These are your protein's residue names, not anyone else's. Find them with:
#
#     python3 scripts/split_receptor.py -i REC.pdb --list
#
# COFACTOR : the organic cofactor to keep.
#            Examples: GDP, GTP, ATP, ADP, NAD, NAP, FAD, FMN, SAM, HEM, COA
#            Set to "" if your protein has no cofactor.
# METAL    : the metal ion to keep.
#            Examples: MG, ZN, MN, CA, FE, CU, NI, K, NA
#            Set to "" if your protein has no metal ion.

COFACTOR="CHANGEME"
COFACTOR_CHARGE=0           # NET FORMAL CHARGE at pH 7, for antechamber -nc.
                            # This one matters: get it wrong and every partial
                            # charge on the molecule is wrong, in a topology
                            # that looks perfectly valid.
                            #   GDP, ADP  -3      GTP, ATP  -4
                            #   NAD+      -1      NADH      -2
                            #   FAD       -2      SAM       +1

METAL="CHANGEME"
METAL_CHARGE=2              # net charge, for reporting only

# The ligand's residue name as it will appear in the topology. Any three
# letters that do not collide with a residue name already in your system.
LIGAND_RESNAME="LIG"
LIGAND_CHARGE=0             # net formal charge. 0 for most neutral drug-like
                            # molecules; -1 for a carboxylate, +1 for an amine
                            # protonated at pH 7.


# --- force field -------------------------------------------------------------
# Numbers as offered by the interactive `gmx pdb2gmx` menu. Run pdb2gmx once by
# hand to confirm them for your GROMACS build - the list order changes between
# versions and this is the one thing that silently picks the wrong force field.
#
#   Common AMBER entries:  1 AMBER03  2 AMBER94  3 AMBER96  4 AMBER99
#                          5 AMBER99SB  6 AMBER99SB-ILDN  7 AMBERGS
#   Water models:          1 TIP3P  2 TIP4P  3 TIP4P-Ew  4 TIP5P
#
# AMBER99SB-ILDN (6) is the usual pairing with GAFF2 ligand parameters and is
# what most protein-ligand literature uses. Change it deliberately, not by
# accident.

FF_CHOICE=6                 # force field menu number
WATER_CHOICE=1              # water model menu number
FF_DIR="amber99sb-ildn.ff"  # must match FF_CHOICE - used for #include paths

# Amber force field files loaded in tleap for the metal ion.
# ions234lm_126_tip3p covers 2+/3+/4+ metals with the 12-6 Lennard-Jones set
# for TIP3P water. Change the water model here if you changed WATER_CHOICE.
ION_FRCMOD="frcmod.ions234lm_126_tip3p"


# --- box and solvent ---------------------------------------------------------

BOX_TYPE="dodecahedron"     # dodecahedron uses ~29% less water than cubic
BOX_DISTANCE=1.0            # nm from solute to box edge. 1.0 is the minimum
                            # that keeps a protein from seeing its own image
                            # with a 1.0 nm cutoff.

WATER_GRO="spc216.gro"      # resolved via GROMACS' own share directory
POSITIVE_ION="NA"
NEGATIVE_ION="CL"
SALT_CONC=0.0               # mol/L of additional salt beyond neutralisation.
                            # Use 0.15 for physiological ionic strength.


# --- simulation --------------------------------------------------------------

MDP_DIR="mdp"               # relative to the repo root
TEMPERATURE=300             # K - must match ref_t in the .mdp files
PRESSURE=1.0                # bar
POSRES_FORCE=1000           # kJ/mol/nm^2 for cofactor and ligand restraints


# --- hardware ----------------------------------------------------------------
# Passed to gmx mdrun for the production run. Leave MDRUN_EXTRA empty to let
# GROMACS decide, which is usually right.
#
#   CPU only          : MDRUN_EXTRA=""
#   One GPU           : MDRUN_EXTRA="-nb gpu -pme gpu -ntmpi 1"
#   Limit CPU threads : MDRUN_EXTRA="-nt 12"

MDRUN_EXTRA=""
CHECKPOINT_MINUTES=15       # -cpt: how often to write md.cpt
MAX_HOURS=0                 # -maxh: wall-clock limit, 0 = no limit


# --- safety ------------------------------------------------------------------
# grompp -maxwarn. Keep this as low as you can bear.
#
# Every warning grompp raises is telling you something real. The usual two on a
# protein-cofactor system are (a) a small non-zero total charge left over from
# the cofactor's fitted partial charges, and (b) an excluded-atom-pair note
# beyond the cutoff. Both are benign. Anything else deserves reading before you
# raise this number.

MAXWARN=2
