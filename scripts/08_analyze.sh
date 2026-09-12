#!/usr/bin/env bash
# Stage 8 - post-processing and standard analyses.
#
# Output: analysis/ containing the centred trajectory, RMSD, RMSF, Rg,
#         hydrogen bonds, energies, and extracted frames.

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

banner "Stage 8/8  -  analysis"

enter_workdir
need_gmx
need_file md.xtc "no production trajectory found - run scripts/07_production.sh"
need_file md.tpr

mkdir -p analysis

# ---------------------------------------------------------------------------
# Periodic boundary treatment. Do this FIRST - every analysis below is wrong
# on a raw trajectory, because molecules jump across the box edge and RMSD
# picks up those jumps as enormous spurious motion.
# ---------------------------------------------------------------------------
if [[ -f analysis/md_center.xtc && "${FORCE:-0}" != "1" ]]; then
    ok "analysis/md_center.xtc exists - skipping trjconv"
else
    step "centring and making molecules whole"
    printf 'Protein\nSystem\n' | \
        gmx trjconv -s md.tpr -f md.xtc -o analysis/md_center.xtc \
                    -center -pbc mol -ur compact \
        || die "trjconv failed. If it cannot find a 'Protein' group, pass an
index file with -n, or use group numbers instead of names."
    ok "analysis/md_center.xtc written"
fi

TRAJ=analysis/md_center.xtc

# ---------------------------------------------------------------------------
step "first and last frames"
printf 'System\n' | gmx trjconv -s md.tpr -f "$TRAJ" \
    -o analysis/frame_first.pdb -dump 0 >/dev/null 2>&1 || true
if [[ -f md.gro ]]; then cp md.gro analysis/frame_last.gro; fi

# ---------------------------------------------------------------------------
step "RMSD  (backbone, fitted on backbone)"
printf 'Backbone\nBackbone\n' | \
    gmx rms -s md.tpr -f "$TRAJ" -o analysis/rmsd.xvg -tu ns >/dev/null 2>&1 \
    && ok "analysis/rmsd.xvg" || warn "RMSD failed - check your group names"

step "RMSD of the ligand relative to the protein"
printf 'Backbone\n%s\n' "$LIGAND_RESNAME" | \
    gmx rms -s md.tpr -f "$TRAJ" -o analysis/rmsd_ligand.xvg -tu ns \
    >/dev/null 2>&1 \
    && ok "analysis/rmsd_ligand.xvg  (fit on protein, measure the ligand -" \
    && echo "        this is the number that tells you whether it stayed bound)" \
    || warn "ligand RMSD failed - you may need an index group for $LIGAND_RESNAME"

# ---------------------------------------------------------------------------
step "RMSF  (per residue)"
printf 'Protein\n' | \
    gmx rmsf -s md.tpr -f "$TRAJ" -o analysis/rmsf.xvg -res >/dev/null 2>&1 \
    && ok "analysis/rmsf.xvg" || warn "RMSF failed"

step "radius of gyration"
printf 'Protein\n' | \
    gmx gyrate -s md.tpr -f "$TRAJ" -o analysis/gyrate.xvg >/dev/null 2>&1 \
    && ok "analysis/gyrate.xvg" || warn "gyrate failed"

# ---------------------------------------------------------------------------
step "hydrogen bonds: protein - ligand"
printf 'Protein\n%s\n' "$LIGAND_RESNAME" | \
    gmx hbond -s md.tpr -f "$TRAJ" -num analysis/hbond_protein_ligand.xvg \
              -tu ns >/dev/null 2>&1 \
    && ok "analysis/hbond_protein_ligand.xvg" \
    || warn "protein-ligand hbond failed (older gmx uses 'gmx hbond', newer 'gmx hbond-legacy')"

if [[ -n "$COFACTOR" ]]; then
    step "hydrogen bonds: protein - $COFACTOR"
    printf 'Protein\n%s\n' "$COFACTOR" | \
        gmx hbond -s md.tpr -f "$TRAJ" \
                  -num "analysis/hbond_protein_${COFACTOR}.xvg" -tu ns \
        >/dev/null 2>&1 \
        && ok "analysis/hbond_protein_${COFACTOR}.xvg" || true
fi

cat <<'EOF'

    A purely hydrophobic ligand has no hydrogen-bond donors or acceptors, so
    zero protein-ligand hydrogen bonds is the correct answer, not a bug. For
    such ligands, judge binding by ligand RMSD, contact counts and buried
    surface area instead.

EOF

# ---------------------------------------------------------------------------
if [[ -n "$METAL" ]]; then
    step "$METAL coordination distances"
    printf 'Protein\n%s\n' "$METAL" | \
        gmx mindist -s md.tpr -f "$TRAJ" -od "analysis/mindist_${METAL}.xvg" \
                    -tu ns >/dev/null 2>&1 \
        && ok "analysis/mindist_${METAL}.xvg  (watch for the ion leaving its site)" \
        || true
fi

# ---------------------------------------------------------------------------
step "energies"
printf 'Potential\nTemperature\nPressure\nDensity\n\n' | \
    gmx energy -f md.edr -o analysis/energy.xvg >/dev/null 2>&1 \
    && ok "analysis/energy.xvg" || warn "energy extraction failed"

# ---------------------------------------------------------------------------
banner "Analysis written to $WORKDIR/analysis"

cat <<EOF
  Plot the .xvg files with xmgrace, or read them with numpy - they are plain
  columns with a header of lines starting @ and #.

  For MM-PBSA / MM-GBSA binding free energies, see docs/WORKFLOW.md. That step
  needs an index file with your protein and ligand groups, and the exact
  gmx_MMPBSA invocation depends on the version you have installed - which is
  why it is documented rather than scripted here.
EOF

note "stage 8 complete: analysis"
