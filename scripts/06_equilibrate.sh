#!/usr/bin/env bash
# Stage 6 - position restraints, then NVT and NPT equilibration.
#
# Output: nvt.gro, npt.gro (+ checkpoints), temperature/pressure/density plots

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

banner "Stage 6/8  -  equilibration (NVT then NPT)"

enter_workdir
need_gmx
need_file em.gro "run scripts/05_minimize.sh first"

# ---------------------------------------------------------------------------
# Position restraints for cofactor and ligand.
#
# pdb2gmx already wrote posre.itp for the protein. The cofactor and ligand
# need their own, or they will drift out of the binding site during the very
# equilibration that is supposed to be holding everything still.
#
# The metal ion is deliberately left unrestrained: it is a single atom whose
# coordination geometry should be allowed to relax with the protein.
# ---------------------------------------------------------------------------
add_posres() {
    local name="$1" itp="${1}_molecule.itp" out="posre_${1}.itp"

    if [[ -f "$out" && "${FORCE:-0}" != "1" ]]; then
        ok "$out exists - skipping"
    else
        need_file "$name.gro"
        step "position restraints for $name (fc = $POSRES_FORCE)"
        # A single-molecule .gro has exactly one selectable group: "System".
        echo 0 | gmx genrestr -f "$name.gro" -o "$out" \
                              -fc "$POSRES_FORCE" "$POSRES_FORCE" "$POSRES_FORCE" \
            || die "gmx genrestr failed for $name"
        n=$(grep -c "$POSRES_FORCE" "$out" || true)
        ok "$n atoms restrained"
    fi

    # Append the #ifdef block to the molecule topology, once.
    if grep -q "posre_${name}.itp" "$itp" 2>/dev/null; then
        ok "$itp already references $out"
    else
        need_file "$itp"
        printf '\n#ifdef POSRES\n#include "%s"\n#endif\n' "$out" >> "$itp"
        ok "appended POSRES block to $itp"
    fi
}

[[ -n "$COFACTOR" ]] && add_posres "$COFACTOR"
add_posres "$LIGAND_RESNAME"
[[ -n "$METAL" ]] && step "$METAL left unrestrained on purpose (single ion)"

# ---------------------------------------------------------------------------
# NVT
# ---------------------------------------------------------------------------
if already_done nvt.gro; then
    ok "NVT already finished"
else
    step "grompp for NVT   (-r is required whenever POSRES is defined)"
    gmx grompp -f "$(mdp nvt.mdp)" -c em.gro -r em.gro -p topol.top \
               -o nvt.tpr -maxwarn "$MAXWARN" 2> grompp_nvt.log \
        || { tail -40 grompp_nvt.log; die "grompp failed for NVT.
If it complains about the tc-grps 'Protein' or 'Non-Protein' not existing,
your system has no protein group under that name - edit tc-grps in
$(mdp nvt.mdp) to match your index groups."; }

    step "mdrun NVT  (500 ps at ${TEMPERATURE} K)"
    gmx mdrun -deffnm nvt ${MDRUN_EXTRA}
    need_file nvt.gro

    echo Temperature | gmx energy -f nvt.edr -o nvt_temperature.xvg 2>/dev/null \
        | grep -iE "^Temperature" || true
fi

# ---------------------------------------------------------------------------
# NPT
# ---------------------------------------------------------------------------
if already_done npt.gro; then
    ok "NPT already finished"
else
    step "grompp for NPT   (-t carries the NVT velocities over)"
    gmx grompp -f "$(mdp npt.mdp)" -c nvt.gro -r nvt.gro -t nvt.cpt \
               -p topol.top -o npt.tpr -maxwarn "$MAXWARN" 2> grompp_npt.log \
        || { tail -40 grompp_npt.log; die "grompp failed for NPT"; }

    step "mdrun NPT  (1 ns at ${TEMPERATURE} K, ${PRESSURE} bar)"
    gmx mdrun -deffnm npt ${MDRUN_EXTRA}
    need_file npt.gro
fi

# ---------------------------------------------------------------------------
# Equilibration report
# ---------------------------------------------------------------------------
step "equilibration summary"
echo
for term in Temperature Pressure Density Volume; do
    out="npt_$(echo "$term" | tr '[:upper:]' '[:lower:]').xvg"
    line=$(echo "$term" | gmx energy -f npt.edr -o "$out" 2>/dev/null \
           | grep -iE "^$term" || true)
    [[ -n "$line" ]] && echo "    $line"
done
echo

cat <<'EOF'
    Check these before committing to a production run:
      Temperature  should sit at your target, fluctuating a few K
      Pressure     fluctuates enormously - only the average means anything
      Density      should be flat by the end, ~1000 kg/m3 for water
      Volume       should stop drifting

    A density still climbing at the end of NPT means the box has not
    equilibrated. Extend NPT rather than starting production.
EOF

note "stage 6 complete: NVT + NPT equilibration"
banner "Stage 6 done  ->  next: scripts/07_production.sh"
