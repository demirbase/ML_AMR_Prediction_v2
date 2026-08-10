#!/usr/bin/env bash
# =============================================================================
# Full $AMR_WORK -> Google Drive backup, in order of what cannot be rebuilt.
# =============================================================================
# Run inside screen/tmux on the login node — a compute node has no internet:
#     screen -S backup
#     bash slurm/rclone_full_backup.sh
#     # detach with Ctrl-A D, reattach with `screen -r backup`
#
# Two hard constraints shape this script:
#
#   1. Google Drive accepts 750 GB per account per 24 h. $AMR_WORK is ~934 GB, so
#      a full copy CANNOT finish in one day. --drive-stop-on-upload-limit makes
#      rclone stop cleanly at the cap instead of thrashing against it; re-run the
#      script the next day and it skips whatever already landed.
#   2. --transfers 4 gets the process OOM-killed on this login node (learned the
#      hard way), so every tier runs at 2.
#
# TIER ORDER IS THE POINT. If the run is cut short — quota, dropped SSH, a node
# going away — what has already been uploaded should still be a usable backup.
# So the irreplaceable ~141 GB goes first and the ~793 GB of regenerable
# intermediates goes last:
#
#   Tier 1  results, models, runs, Zenodo tarballs   outputs, nothing rebuilds them
#   Tier 2  containers, data/external                the pinned software + CARD snapshot
#   Tier 3  data/processed/*/lineage                 PopPUNK clusters, expensive to redo
#   Tier 4  backup/                                  older on-HPC backup, kept as-is
#   Tier 5  data/raw                                 17.7k genomes; BV-BRC drifts, so back up
#   Tier 6  data/processed/*/{unitig_all,<ab>}       matrices + stores: deterministic
#                                                    from tiers 2+5, just slow
#
# Re-running is safe and is the intended way to finish: rclone copy compares size
# and modtime and skips what matches.
# =============================================================================
set -uo pipefail            # deliberately NOT -e: a tier hitting the daily cap
                            # must not abort the tiers that follow next run

: "${AMR_WORK:?AMR_WORK is not set — export it first}"
DEST="${DEST:-gdrive:AMR_FULL_20260805}"
LOG="${LOG:-$AMR_WORK/rclone_full_backup.log}"

RCLONE_OPTS=(
  --transfers 2 --checkers 4
  --drive-stop-on-upload-limit          # stop at the 750 GB/day cap, do not thrash
  --retries 5 --low-level-retries 20
  --stats 5m --stats-one-line
  --log-file "$LOG" --log-level INFO
)

copy_tier () {                          # copy_tier <label> <relative path>
  local label="$1" rel="$2"
  local src="$AMR_WORK/$rel"
  [[ -e "$src" ]] || { echo "  skip  $label ($rel not present)"; return 0; }
  local size; size=$(du -sh "$src" 2>/dev/null | cut -f1)
  echo "  ==>   $label  [$rel, $size]  $(date '+%H:%M:%S')"
  rclone copy "$src" "$DEST/$rel" "${RCLONE_OPTS[@]}"
  local rc=$?
  case $rc in
    0) echo "  done  $label" ;;
    7) echo "  CAP   $label — daily upload limit reached; re-run tomorrow to continue"
       echo "$(date -Iseconds) HIT DAILY LIMIT during $label" >> "$LOG"
       exit 7 ;;
    *) echo "  FAIL  $label (rclone exit $rc) — continuing to the next tier" ;;
  esac
}

echo "=============================================================="
echo " FULL BACKUP  $AMR_WORK  ->  $DEST"
echo " started $(date -Iseconds)   ·   log: $LOG"
echo "=============================================================="

# --- Tier 1: the delivered artefacts ----------------------------------------
copy_tier "results (KB, tables, figures, per-model outputs)" "results"
copy_tier "models"                                            "models"
copy_tier "runs (provenance)"                                 "runs"
copy_tier "logs"                                              "logs"
for f in amrk-db_v0.7.1.tar.gz amrk-db_v0.7.1_config-docs.tar.gz; do
  [[ -f "$AMR_WORK/$f" ]] && {
    echo "  ==>   Zenodo tarball $f"
    rclone copy "$AMR_WORK/$f" "$DEST/" "${RCLONE_OPTS[@]}"
  }
done

# --- Tier 2: the environment the results depend on --------------------------
copy_tier "containers (.sif — bit-exact software)"            "containers"
copy_tier "data/external (CARD snapshot etc.)"                "data/external"

# --- Tier 3: population structure -------------------------------------------
for org in "$AMR_WORK"/data/processed/*/lineage; do
  [[ -d "$org" ]] || continue
  rel="${org#$AMR_WORK/}"
  copy_tier "lineage $(basename "$(dirname "$org")")"         "$rel"
done

# --- Tier 4 / 5: older backup + raw genomes ---------------------------------
copy_tier "backup (older on-HPC copy)"                        "backup"
copy_tier "data/raw (17.7k genomes)"                          "data/raw"

# --- Tier 6: regenerable intermediates (the ~793 GB tail) -------------------
echo
echo "  --- tier 6: regenerable intermediates; expect the daily cap here ---"
for d in "$AMR_WORK"/data/processed/*/*; do
  [[ -d "$d" ]] || continue
  [[ "$(basename "$d")" == "lineage" ]] && continue      # already done in tier 3
  rel="${d#$AMR_WORK/}"
  copy_tier "$(basename "$(dirname "$d")")/$(basename "$d")" "$rel"
done

echo "=============================================================="
echo " finished $(date -Iseconds)"
echo " verify with:  rclone check \"$AMR_WORK/results\" \"$DEST/results\" --size-only"
echo "=============================================================="
