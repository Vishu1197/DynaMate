#!/usr/bin/env bash
# Stage 7 - unrestrained production MD.
#
# Safe to re-run: if md.cpt exists this resumes from the checkpoint rather than
# starting over. That is the behaviour you want after a power cut, and the
# behaviour you must NOT get wrong - regenerating md.tpr throws the run away.
#
# Output: md.xtc, md.edr, md.log, md.gro

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

banner "Stage 7/8  -  production molecular dynamics"

enter_workdir
need_gmx
need_file npt.gro "run scripts/06_equilibrate.sh first"

# ---------------------------------------------------------------------------
# Build the tpr once and only once.
# ---------------------------------------------------------------------------
if [[ -f md.tpr ]]; then
    ok "md.tpr already exists - reusing it (this is what makes resume work)"
else
    step "grompp for production"
    warn "POSRES is deliberately NOT defined here. Check $(mdp md.mdp)"
    warn "has 'define =' with nothing after it - a stray -DPOSRES means you"
    warn "simulate a frozen protein for the entire run."
    gmx grompp -f "$(mdp md.mdp)" -c npt.gro -t npt.cpt -p topol.top \
               -o md.tpr -maxwarn "$MAXWARN" 2> grompp_md.log \
        || { tail -40 grompp_md.log; die "grompp failed for production"; }

    nsteps=$(grep -oP 'nsteps\s*=\s*\K[0-9]+' "$(mdp md.mdp)" | head -1)
    dt=$(grep -oP '^\s*dt\s*=\s*\K[0-9.]+' "$(mdp md.mdp)" | head -1)
    if [[ -n "$nsteps" && -n "$dt" ]]; then
        ns=$(awk "BEGIN{printf \"%.1f\", $nsteps * $dt / 1000}")
        ok "production length: $nsteps steps x $dt ps = $ns ns"
    fi
fi

# ---------------------------------------------------------------------------
# Run, resuming if a checkpoint is present.
# ---------------------------------------------------------------------------
run_args=(-s md.tpr -deffnm md -cpo md.cpt -cpt "$CHECKPOINT_MINUTES"
          -maxh "$MAX_HOURS")
[[ -n "$MDRUN_EXTRA" ]] && read -r -a extra <<< "$MDRUN_EXTRA" && run_args+=("${extra[@]}")

if [[ -f md.cpt ]]; then
    step "md.cpt found - RESUMING from the checkpoint"
    ok "no work is lost; the run continues toward the original endpoint"
    run_args+=(-cpi md.cpt)
else
    step "starting production from npt.gro"
fi

echo
echo "    gmx mdrun ${run_args[*]}"
echo

gmx mdrun "${run_args[@]}"

# ---------------------------------------------------------------------------
if [[ -f md.gro ]]; then
    ok "production complete - md.gro written"
    note "stage 7 complete: production MD finished"
    banner "Stage 7 done  ->  next: scripts/08_analyze.sh"
else
    warn "md.gro not present, so the run has not reached its endpoint."
    warn "If it stopped early (wall-clock limit, interruption, crash), simply"
    warn "run this script again - it will pick up from md.cpt."
    note "stage 7 partial: production interrupted"
fi
