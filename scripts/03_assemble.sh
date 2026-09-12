#!/usr/bin/env bash
# Stage 3 - protein topology, then assemble every component into one system.
#
# Output: topol.top (patched), complex.gro

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

banner "Stage 3/8  -  protein topology and system assembly"

enter_workdir
need_gmx
need_file protein.pdb "run scripts/01_prepare.sh first"

# ---------------------------------------------------------------------------
# Protein topology
# ---------------------------------------------------------------------------
if [[ -f protein_processed.gro && "${FORCE:-0}" != "1" ]]; then
    ok "protein_processed.gro exists - skipping pdb2gmx (FORCE=1 to redo)"
else
    step "pdb2gmx  (force field $FF_CHOICE, water model $WATER_CHOICE)"
    warn "confirm those menu numbers match your GROMACS build - the list"
    warn "order changes between versions and a wrong pick is silent."
    printf '%s\n%s\n' "$FF_CHOICE" "$WATER_CHOICE" | \
        gmx pdb2gmx -f protein.pdb -o protein_processed.gro -p topol.top \
                    -i posre.itp -ignh \
        || die "pdb2gmx failed. Common causes: non-standard residues still in \
protein.pdb, missing backbone atoms, or a chain break needing -ter handling."
fi

need_file topol.top
need_file protein_processed.gro

nat=$(sed -n '2p' protein_processed.gro | tr -d ' ')
ok "protein: $nat atoms"

# ---------------------------------------------------------------------------
# Assemble coordinates. ORDER MATTERS - it must match [ molecules ] exactly.
# ---------------------------------------------------------------------------
step "merging coordinates"
merge_inputs=(protein_processed.gro)
[[ -n "$COFACTOR" ]] && { need_file "$COFACTOR.gro"; merge_inputs+=("$COFACTOR.gro"); }
[[ -n "$METAL"    ]] && { need_file "$METAL.gro";    merge_inputs+=("$METAL.gro"); }
need_file "$LIGAND_RESNAME.gro"
merge_inputs+=("$LIGAND_RESNAME.gro")

python3 "$(helper merge_gro.py)" -o complex.gro \
        --title "$SYSTEM_NAME" \
        --box-from protein_processed.gro \
        "${merge_inputs[@]}"

need_file complex.gro

# ---------------------------------------------------------------------------
# Patch topol.top
# ---------------------------------------------------------------------------
step "patching topol.top"

at_args=(); mol_args=(); molecules_args=()

protein_mol=$(awk '/^\[ *moleculetype *\]/{f=1;next}
                   f && !/^;/ && NF {print $1; exit}' topol.top)
protein_mol="${protein_mol:-Protein_chain_A}"
molecules_args+=(--molecules "${protein_mol}:1")

if [[ -n "$COFACTOR" ]]; then
    at_args+=(--atomtypes "${COFACTOR}_atomtypes.itp")
    mol_args+=(--molecule "${COFACTOR}_molecule.itp")
    molecules_args+=(--molecules "${COFACTOR}:1")
fi
if [[ -n "$METAL" ]]; then
    # Atom types only. The molecule itself comes from the force field's ions.itp.
    [[ -f "${METAL}_atomtypes.itp" ]] && at_args+=(--atomtypes "${METAL}_atomtypes.itp")
    molecules_args+=(--molecules "${METAL}:1")
fi
at_args+=(--atomtypes "${LIGAND_RESNAME}_atomtypes.itp")
mol_args+=(--molecule "${LIGAND_RESNAME}_molecule.itp")
molecules_args+=(--molecules "${LIGAND_RESNAME}:1")

python3 "$(helper build_topology.py)" -p topol.top \
        "${at_args[@]}" "${mol_args[@]}" "${molecules_args[@]}"

# ---------------------------------------------------------------------------
# Validate before solvating. Catching a topology error now costs seconds;
# catching it after solvation costs a rebuild.
# ---------------------------------------------------------------------------
step "validating the unsolvated system with grompp"
if gmx grompp -f "$(mdp em.mdp)" -c complex.gro -p topol.top \
              -o validate.tpr -maxwarn "$MAXWARN" 2> grompp_validate.log; then
    ok "topology and coordinates agree"
    grep -iE "System has non-zero total charge" grompp_validate.log || true
    rm -f validate.tpr
else
    tail -40 grompp_validate.log
    die "grompp rejected the unsolvated system. The usual causes:
   - [ molecules ] order does not match the .gro atom order
   - a duplicate atom type (you included both NAME.top and NAME_atomtypes.itp)
   - the metal's moleculetype included twice (force field ions.itp + your own)
   Full output is in $WORKDIR/grompp_validate.log"
fi

note "stage 3 complete: complex.gro assembled and validated"
banner "Stage 3 done  ->  next: scripts/04_solvate.sh"
