#!/usr/bin/env bash
# Stage 2 - generate GAFF2 parameters for the cofactor and ligand, and Amber
#           ion parameters for the metal, then convert all of them to GROMACS.
#
# This is the slowest and most fragile stage. antechamber's AM1-BCC charge
# fitting on a cofactor the size of GDP can take several minutes.
#
# Output: <NAME>_atomtypes.itp, <NAME>_molecule.itp, <NAME>.gro  for each component

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

banner "Stage 2/8  -  parameterise cofactor, metal and ligand"

enter_workdir
need_ambertools
need_cmd tleap "AmberTools is not fully on PATH"

AMB2GRO=""
for cand in amb2gro_top_gro.py acpype; do
    if command -v "$cand" >/dev/null 2>&1; then AMB2GRO="$cand"; break; fi
done
[[ -n "$AMB2GRO" ]] || die "neither amb2gro_top_gro.py nor acpype found. \
Both ship with AmberTools; install one before continuing."

# ---------------------------------------------------------------------------
# An organic molecule: antechamber -> parmchk2 -> tleap -> GROMACS
# $1 = residue/base name, $2 = input pdb, $3 = net charge
# ---------------------------------------------------------------------------
parameterize_organic() {
    local name="$1" pdb="$2" charge="$3"

    if [[ -f "${name}_molecule.itp" && "${FORCE:-0}" != "1" ]]; then
        ok "${name}_molecule.itp exists - skipping (FORCE=1 to redo)"
        return
    fi

    need_file "$pdb"
    step "$name: AM1-BCC charges (net charge $charge) - this can take minutes"
    antechamber -i "$pdb" -fi pdb -o "$name.mol2" -fo mol2 \
                -at gaff2 -c bcc -nc "$charge" -rn "$name" -pf y \
        || die "antechamber failed on $pdb. Check that the molecule has \
hydrogens, sensible geometry, and that the net charge $charge is right."

    step "$name: filling parameter gaps with parmchk2"
    parmchk2 -i "$name.mol2" -f mol2 -o "$name.frcmod" -s 2

    step "$name: building Amber topology with tleap"
    cat > "leap_${name}.in" <<EOF
source leaprc.gaff2
loadamberparams $name.frcmod
$name = loadmol2 $name.mol2
check $name
charge $name
saveamberparm $name $name.prmtop $name.inpcrd
quit
EOF
    tleap -f "leap_${name}.in" > "leap_${name}.log" 2>&1 \
        || { tail -30 "leap_${name}.log"; die "tleap failed for $name"; }

    if grep -qiE "fatal|could not|failed" "leap_${name}.log"; then
        warn "tleap reported problems for $name - read leap_${name}.log"
        grep -iE "fatal|could not|failed|warning" "leap_${name}.log" | head -10
    fi
    need_file "$name.prmtop" "tleap did not produce a topology"

    step "$name: converting to GROMACS format"
    convert_to_gromacs "$name"

    step "$name: extracting includable .itp files"
    python3 "$(helper extract_itp.py)" -i "$name.top" -n "$name"

    ok "$name parameterised"
}

convert_to_gromacs() {
    local name="$1"
    if [[ "$AMB2GRO" == "amb2gro_top_gro.py" ]]; then
        amb2gro_top_gro.py -p "$name.prmtop" -c "$name.inpcrd" \
                           -t "$name.top" -g "$name.gro" \
                           -b "${name}_converted.pdb" \
            || die "amb2gro_top_gro.py failed for $name"
    else
        acpype -p "$name.prmtop" -x "$name.inpcrd" -b "$name" \
            || die "acpype failed for $name"
        cp "${name}.amb2gmx/${name}_GMX.top" "$name.top"
        cp "${name}.amb2gmx/${name}_GMX.gro" "$name.gro"
    fi
    need_file "$name.top"
    need_file "$name.gro"
}

# ---------------------------------------------------------------------------
# Cofactor
# ---------------------------------------------------------------------------
if [[ -n "$COFACTOR" ]]; then
    parameterize_organic "$COFACTOR" "$COFACTOR.pdb" "$COFACTOR_CHARGE"
else
    step "no COFACTOR set - skipping"
fi

# ---------------------------------------------------------------------------
# Ligand
# ---------------------------------------------------------------------------
parameterize_organic "$LIGAND_RESNAME" "$LIGAND_PDB" "$LIGAND_CHARGE"

# ---------------------------------------------------------------------------
# Metal ion - a monatomic ion needs no GAFF treatment, only the Amber ion
# parameter set. And its [ moleculetype ] is ALREADY in the force field's
# ions.itp, so we take the atom types only. Including both is a duplicate
# definition error that is very hard to read.
# ---------------------------------------------------------------------------
if [[ -n "$METAL" ]]; then
    if [[ -f "${METAL}_atomtypes.itp" && "${FORCE:-0}" != "1" ]]; then
        ok "${METAL}_atomtypes.itp exists - skipping"
    else
        need_file "$METAL.pdb"
        step "$METAL: Amber ion parameters ($ION_FRCMOD)"
        cat > "leap_${METAL}.in" <<EOF
source leaprc.gaff2
loadamberparams $ION_FRCMOD
loadOff atomic_ions.lib
$METAL = loadpdb $METAL.pdb
check $METAL
charge $METAL
saveamberparm $METAL $METAL.prmtop $METAL.inpcrd
quit
EOF
        tleap -f "leap_${METAL}.in" > "leap_${METAL}.log" 2>&1 \
            || { tail -30 "leap_${METAL}.log"; die "tleap failed for $METAL"; }

        charge_line=$(grep -i "Total unperturbed charge" "leap_${METAL}.log" | head -1 || true)
        [[ -n "$charge_line" ]] && ok "$METAL $charge_line"

        convert_to_gromacs "$METAL"
        python3 "$(helper extract_itp.py)" -i "$METAL.top" -n "$METAL" \
                --atomtypes-only
        ok "$METAL parameterised (atom types only - the force field already \
defines the ion itself)"
    fi
else
    step "no METAL set - skipping"
fi

note "stage 2 complete: parameterised ${COFACTOR:-none} / ${METAL:-none} / $LIGAND_RESNAME"
banner "Stage 2 done  ->  next: scripts/03_assemble.sh"
