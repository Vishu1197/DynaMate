#!/usr/bin/env bash
# Stage 4 - define the box, add water, neutralise with counter-ions.
#
# Output: solvated_ions.gro, topol.top updated with SOL and ion counts

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

banner "Stage 4/8  -  solvation and counter-ions"

enter_workdir
need_gmx
need_file complex.gro "run scripts/03_assemble.sh first"
need_file topol.top

step "defining a $BOX_TYPE box, $BOX_DISTANCE nm clearance"
gmx editconf -f complex.gro -o boxed.gro -c -d "$BOX_DISTANCE" -bt "$BOX_TYPE"

step "adding water"
gmx solvate -cp boxed.gro -cs "$WATER_GRO" -p topol.top -o solvated.gro \
    || die "gmx solvate failed. If it cannot find $WATER_GRO, source GMXRC so \
GROMACS can locate its own share directory."

nsol=$(awk '/^\[ *molecules *\]/{f=1;next} f && $1=="SOL" {print $2}' topol.top | tail -1)
ok "${nsol:-?} water molecules added"

step "building a tpr for genion"
gmx grompp -f "$(mdp ions.mdp)" -c solvated.gro -p topol.top \
           -o ions.tpr -maxwarn "$MAXWARN" 2> grompp_ions.log \
    || { tail -30 grompp_ions.log; die "grompp failed before ion placement"; }

charge=$(grep -oP 'System has non-zero total charge: \K[-0-9.]+' grompp_ions.log | tail -1 || true)
if [[ -n "$charge" ]]; then
    ok "system charge before neutralisation: $charge"
else
    ok "system is already neutral"
fi

step "replacing water with counter-ions"
genion_args=(-s ions.tpr -o solvated_ions.gro -p topol.top
             -pname "$POSITIVE_ION" -nname "$NEGATIVE_ION" -neutral)
if awk "BEGIN{exit !($SALT_CONC > 0)}"; then
    genion_args+=(-conc "$SALT_CONC")
    step "  plus $SALT_CONC mol/L background salt"
fi

# genion asks which group to take solvent from. SOL is what you want -
# answering with anything else replaces protein atoms with ions.
echo SOL | gmx genion "${genion_args[@]}" \
    || die "gmx genion failed"

need_file solvated_ions.gro

step "final composition"
awk '/^\[ *molecules *\]/{f=1} f' topol.top | sed 's/^/    /'

warn "read that table. Every component you expect should be there, exactly once."
warn "A missing cofactor here means it was never in [ molecules ] and grompp"
warn "will keep working while simulating a system you did not intend."

note "stage 4 complete: solvated and neutralised"
banner "Stage 4 done  ->  next: scripts/05_minimize.sh"
