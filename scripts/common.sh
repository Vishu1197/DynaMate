#!/usr/bin/env bash
# Shared helpers for the DynaMate stage scripts. Sourced, not run.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ ! -f "$REPO_ROOT/config.sh" ]]; then
    echo "config.sh not found in $REPO_ROOT" >&2
    exit 1
fi
# shellcheck disable=SC1091
source "$REPO_ROOT/config.sh"

# ---------------------------------------------------------------------------
# Refuse to run on an unedited config.
#
# DynaMate ships with CHANGEME sentinels rather than a working example, so that
# nobody can clone it and accidentally run someone else's system. If you see
# this message, open config.sh - that is the whole setup.
# ---------------------------------------------------------------------------
_unset=()
for _v in WORKDIR COFACTOR METAL; do
    if [[ "${!_v-}" == "CHANGEME" ]]; then _unset+=("$_v"); fi
done
if (( ${#_unset[@]} )); then
    printf '\n\033[1;31mconfig.sh has not been filled in.\033[0m\n\n' >&2
    printf '  Still set to CHANGEME: %s\n\n' "${_unset[*]}" >&2
    cat >&2 <<'EOF'
  Open config.sh and set them for YOUR system.

    WORKDIR   an absolute path where everything will be created
    COFACTOR  the cofactor's residue name in your PDB, or "" if there is none
    METAL     the metal ion's residue name in your PDB, or "" if there is none

  To find the residue names in your own structure:

    python3 scripts/split_receptor.py -i /path/to/REC.pdb --list

  It prints every residue in the file with its atom count. Use those names.

EOF
    exit 1
fi
unset _unset _v

# ---------------------------------------------------------------- output ----

C_HEAD=$'\033[1;36m'; C_OK=$'\033[1;32m'; C_WARN=$'\033[1;33m'
C_ERR=$'\033[1;31m';  C_OFF=$'\033[0m'
if [[ ! -t 1 ]]; then C_HEAD=""; C_OK=""; C_WARN=""; C_ERR=""; C_OFF=""; fi

banner() { printf '%s\n%s  %s\n%s\n' \
    "${C_HEAD}============================================================================" \
    "" "$*" "============================================================================${C_OFF}"; }
step()   { printf '%s>>%s %s\n' "$C_HEAD" "$C_OFF" "$*"; }
ok()     { printf '%s ok%s %s\n' "$C_OK" "$C_OFF" "$*"; }
warn()   { printf '%s  !%s %s\n' "$C_WARN" "$C_OFF" "$*" >&2; }
die()    { printf '%serror%s %s\n' "$C_ERR" "$C_OFF" "$*" >&2; exit 1; }

# ----------------------------------------------------------------- checks ----

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "'$1' not found in PATH. $2"
}

need_gmx() {
    need_cmd gmx "Source your GROMACS install, e.g. 'source /usr/local/gromacs/bin/GMXRC'"
}

need_ambertools() {
    command -v antechamber >/dev/null 2>&1 || die \
        "antechamber not found. Activate your AmberTools environment first, e.g. 'conda activate ambertools'"
}

need_file() {
    [[ -f "$1" ]] || die "expected file not found: $1${2:+  ($2)}"
}

# Enter the working directory, creating it if necessary.
enter_workdir() {
    mkdir -p "$WORKDIR"
    cd "$WORKDIR"
}

# Absolute path to an mdp file in the repo.
mdp() { echo "$REPO_ROOT/$MDP_DIR/$1"; }

# Absolute path to a helper script in the repo.
helper() { echo "$REPO_ROOT/scripts/$1"; }

# Refuse to clobber a finished stage unless FORCE=1.
already_done() {
    local marker="$1"
    if [[ -f "$marker" && "${FORCE:-0}" != "1" ]]; then
        ok "$marker already exists - skipping. Set FORCE=1 to redo this stage."
        return 0
    fi
    return 1
}

# Record a one-line note in the project log.
note() {
    mkdir -p "$WORKDIR"
    printf '%s  %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >> "$WORKDIR/dynamate.log"
}
