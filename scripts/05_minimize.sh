#!/usr/bin/env bash
# Stage 5 - energy minimisation.
#
# Output: em.gro, em.edr, potential.xvg

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

banner "Stage 5/8  -  energy minimisation"

enter_workdir
need_gmx
need_file solvated_ions.gro "run scripts/04_solvate.sh first"

if already_done em.gro; then
    banner "Stage 5 done  ->  next: scripts/06_equilibrate.sh"; exit 0
fi

step "grompp"
gmx grompp -f "$(mdp em.mdp)" -c solvated_ions.gro -p topol.top \
           -o em.tpr -maxwarn "$MAXWARN" 2> grompp_em.log \
    || { tail -40 grompp_em.log; die "grompp failed for energy minimisation"; }

step "mdrun  (steepest descent, this is quick)"
gmx mdrun -v -deffnm em

need_file em.gro

step "checking convergence"
echo Potential | gmx energy -f em.edr -o potential.xvg >/dev/null 2>&1 || true

fmax=$(grep -i "Maximum force" em.log | tail -1 || true)
pot=$(grep -i "Potential Energy" em.log | tail -1 || true)
steps=$(grep -iE "converged to Fmax|did not converge" em.log | tail -1 || true)

echo
[[ -n "$steps" ]] && echo "    $steps"
[[ -n "$pot"   ]] && echo "    $pot"
[[ -n "$fmax"  ]] && echo "    $fmax"
echo

if grep -qi "did not converge" em.log; then
    warn "minimisation did NOT reach the force tolerance."
    warn "That is not automatically fatal, but a large Fmax usually means a"
    warn "bad contact - most often the docked ligand overlapping the protein,"
    warn "or a cofactor placed from a different coordinate frame."
    warn "Look at em.gro before continuing."
else
    ok "converged"
fi

if [[ -n "$pot" ]]; then
    val=$(echo "$pot" | grep -oE '[-0-9.e+]+' | tail -1)
    if awk "BEGIN{exit !($val > 0)}" 2>/dev/null; then
        warn "potential energy is positive - the system is still strained."
    fi
fi

note "stage 5 complete: energy minimisation"
banner "Stage 5 done  ->  next: scripts/06_equilibrate.sh"
