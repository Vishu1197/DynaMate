#!/usr/bin/env bash
# Stage 1 - split the receptor into protein, cofactor and metal ion.
#
# Input : $RECEPTOR_PDB  (protein + cofactor + metal, from the crystal structure)
#         $LIGAND_PDB    (the docked ligand alone, H-added and minimised)
# Output: protein.pdb, <COFACTOR>.pdb, <METAL>.pdb

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

banner "Stage 1/8  -  prepare input structures"

enter_workdir
need_file "$RECEPTOR_PDB" "put your receptor in $WORKDIR"
need_file "$LIGAND_PDB"   "put your ligand in $WORKDIR"

step "inventory of $RECEPTOR_PDB"
python3 "$(helper split_receptor.py)" -i "$RECEPTOR_PDB" --list

extract_args=()
[[ -n "$COFACTOR" ]] && extract_args+=(--extract "$COFACTOR")
[[ -n "$METAL"    ]] && extract_args+=(--extract "$METAL")

step "splitting"
python3 "$(helper split_receptor.py)" -i "$RECEPTOR_PDB" "${extract_args[@]}"

need_file protein.pdb
[[ -n "$COFACTOR" ]] && need_file "$COFACTOR.pdb"
[[ -n "$METAL"    ]] && need_file "$METAL.pdb"

step "checking the ligand has hydrogens"
h_count=$(awk '/^(ATOM|HETATM)/ {e=substr($0,77,2); gsub(/ /,"",e);
                if (e=="H" || substr($0,14,1)=="H") n++} END {print n+0}' \
          "$LIGAND_PDB")
if [[ "$h_count" -eq 0 ]]; then
    warn "$LIGAND_PDB appears to contain no hydrogen atoms."
    warn "antechamber needs an all-atom ligand. Add hydrogens and minimise it"
    warn "first (Chimera: Tools > Structure Editing > AddH, then Minimize)."
else
    ok "$h_count hydrogens found in the ligand"
fi

note "stage 1 complete: split $RECEPTOR_PDB"
banner "Stage 1 done  ->  next: scripts/02_parameterize.sh"
