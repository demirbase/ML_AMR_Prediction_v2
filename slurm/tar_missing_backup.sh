#!/usr/bin/env bash
# =============================================================================
# Repair pass for the full backup: the two tiers it lost to SIGKILL.
# =============================================================================
# What happened on 2026-08-11..13: rclone_full_backup.sh ships file-by-file. On
# the login node the long, many-small-file transfers were killed (exit 137):
#
#   results     18,621 of 20,365 files never landed  (17,814 of them
#               {org}/global_exploration — exactly the "many tiny files" shape)
#   data/raw    17,742 genomes, nothing landed at all
#
# The script maps every exit code other than 0 and 7 to "continue to the next
# tier" and then prints `finished` unconditionally, so the hole stayed invisible
# for five days. This script fixes the transfer shape, not the symptom: stop
# shipping 36k small files, tar each subtree into one archive and ship a handful
# of large files. rclone is fast and stable on those, and each upload is short
# enough that nothing gets killed halfway.
#
# Two nodes are needed, because of two hard facts about this cluster:
#   tar on the login node dies with SIGXCPU  -> tar runs on a debug node
#   compute nodes have no internet           -> rclone runs on the login node
#
# The phases therefore alternate, and staging disk is reused between groups so
# peak usage stays at one group instead of both:
#
#   # 1. debug node
#   srun -p debug -N1 -c8 --time=01:00:00 --pty bash
#   bash slurm/tar_missing_backup.sh tar results
#   exit
#   # 2. login node, inside screen (Ctrl-A D to detach)
#   screen -S tarup
#   bash slurm/tar_missing_backup.sh upload results
#   bash slurm/tar_missing_backup.sh clean  results    # only after upload verified
#   # 3. same three steps for the genomes
#   bash slurm/tar_missing_backup.sh tar raw           # on the debug node
#   bash slurm/tar_missing_backup.sh upload raw        # on the login node
#   bash slurm/tar_missing_backup.sh clean  raw
#
# Everything lands in ONE Drive folder ($DEST, default the existing
# AMR_FULL_20260805) under tar/. Those tarballs are the authoritative copy of
# `results` and `data/raw`; the file-level results/ tree already up there is the
# partial one, left in place because it costs nothing.
#
# MANIFEST.tsv records, per archive: source path, source file count, source
# bytes, archive bytes, MD5. That is what answers "is this archive complete"
# later, for Zenodo or for the defence — a tarball alone cannot answer it.
# =============================================================================
set -uo pipefail

: "${AMR_WORK:?AMR_WORK is not set — export it first}"
STAGE="${STAGE:-$AMR_WORK/backup_tars}"
DEST="${DEST:-gdrive:AMR_FULL_20260805/tar}"
MANIFEST="$STAGE/MANIFEST.tsv"
LOG="${LOG:-$AMR_WORK/tar_missing_backup.log}"

PHASE="${1:-}"; GROUP="${2:-}"
case "$PHASE/$GROUP" in
  tar/results|tar/raw|upload/results|upload/raw|clean/results|clean/raw) ;;
  *) echo "usage: $0 {tar|upload|clean} {results|raw}" >&2; exit 2 ;;
esac

case "$GROUP" in
  results) ROOT="$AMR_WORK/results" ;;
  raw)     ROOT="$AMR_WORK/data/raw" ;;
esac
[[ -d "$ROOT" ]] || { echo "no such source: $ROOT" >&2; exit 2; }

# Subtrees = immediate subdirectories, so one dead upload costs one subtree
# rather than the whole group. A flat directory becomes a single archive.
mapfile -t SUBS < <(find "$ROOT" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
FLAT=0
if [[ ${#SUBS[@]} -eq 0 ]]; then FLAT=1; SUBS=("$(basename "$ROOT")"); fi

# gzip -1: the trees are a mix of text (csv/json/txt, compresses well) and
# already-compressed npz/png/pdf (does not). Level 1 keeps the CPU cost near
# zero on the parts that will not shrink anyway.
if command -v pigz >/dev/null 2>&1; then
  COMPRESS=(pigz -1 -p "${SLURM_CPUS_PER_TASK:-8}"); EXT="tar.gz"
else
  COMPRESS=(gzip -1); EXT="tar.gz"
fi

log () { echo "$(date '+%F %T')  $*" | tee -a "$LOG"; }

# --- tar ---------------------------------------------------------------------
do_tar () {
  mkdir -p "$STAGE"
  local avail_kb; avail_kb=$(df -Pk "$STAGE" | awk 'NR==2{print $4}')
  log "stage $STAGE — $((avail_kb/1024/1024)) GiB free; source $ROOT is $(du -sh "$ROOT" 2>/dev/null | cut -f1)"
  [[ -f "$MANIFEST" ]] || printf 'archive\tsource\tn_files\tsrc_bytes\ttar_bytes\tmd5\n' > "$MANIFEST"

  local sub src out n sb tb md5 rc
  for sub in "${SUBS[@]}"; do
    if [[ $FLAT -eq 1 ]]; then src="$ROOT"; else src="$ROOT/$sub"; fi
    out="$STAGE/${GROUP}__${sub}.${EXT}"
    # Skip only what is genuinely current. "Already built" was the wrong test: after
    # results/tables and results/figures were regenerated on this machine, a re-run
    # skipped both and the backup would have silently kept the stale archives. Compare
    # against the source instead -- one file newer than the archive means rebuild.
    if [[ -s "$out" ]] && grep -q "^$(basename "$out")	" "$MANIFEST"; then
      if [[ -z "$(find "$src" -type f -newer "$out" -print -quit 2>/dev/null)" ]]; then
        log "skip  $(basename "$out") (archive newer than every source file)"; continue
      fi
      log "stale $(basename "$out") — source changed since it was built; rebuilding"
      # Drop the old row first: appending would leave two rows for one archive and
      # `rclone check` would then compare against whichever the reader hit first.
      grep -v "^$(basename "$out")	" "$MANIFEST" > "$MANIFEST.tmp" && mv "$MANIFEST.tmp" "$MANIFEST"
      rm -f "$STAGE/.verified_$GROUP"     # a rebuilt archive is not a verified archive
    fi
    n=$(find "$src" -type f | wc -l)
    sb=$(du -sb "$src" | cut -f1)
    log "tar   $(basename "$out")  <- $sub  ($n files, $((sb/1024/1024)) MiB)"
    tar -C "$(dirname "$src")" -cf - "$(basename "$src")" | "${COMPRESS[@]}" > "$out"
    rc=${PIPESTATUS[0]}${PIPESTATUS[1]}
    if [[ "$rc" != "00" ]]; then
      log "FAIL  $(basename "$out") (tar/compress exit $rc) — removing partial archive"
      rm -f "$out"; continue
    fi
    tb=$(stat -c%s "$out")
    md5=$(md5sum "$out" | cut -d' ' -f1)
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(basename "$out")" "${src#$AMR_WORK/}" "$n" "$sb" "$tb" "$md5" >> "$MANIFEST"
    log "done  $(basename "$out")  $((tb/1024/1024)) MiB  md5 $md5"
  done
  log "tar phase done for '$GROUP'; manifest: $MANIFEST"
}

# --- upload ------------------------------------------------------------------
# One file per rclone invocation: a short-lived process cannot be killed for
# running long, and a death costs one archive instead of the whole group.
# --buffer-size 0 --use-mmap keep the footprint small; that footprint is what
# got the original run SIGKILLed.
do_upload () {
  [[ -d "$STAGE" ]] || { echo "nothing staged in $STAGE — run the tar phase first" >&2; exit 2; }
  local opts=(--transfers 1 --checkers 2 --buffer-size 0 --use-mmap
              --drive-chunk-size 32M --drive-stop-on-upload-limit
              --retries 5 --low-level-retries 20
              --stats 2m --stats-one-line --log-file "$LOG" --log-level INFO)
  local failed=() f rc
  for f in "$STAGE/${GROUP}__"*.${EXT}; do
    [[ -e "$f" ]] || { log "no archives for '$GROUP' in $STAGE"; return 1; }
    log "up    $(basename "$f")  $(( $(stat -c%s "$f") /1024/1024 )) MiB"
    rclone copy "$f" "$DEST/" "${opts[@]}"; rc=$?
    case $rc in
      0) log "done  $(basename "$f")" ;;   # rclone skips size+modtime matches; the md5 check below is the real gate
      7) log "CAP   daily 750 GB limit reached at $(basename "$f") — re-run this same command tomorrow"
         return 7 ;;
      *) log "FAIL  $(basename "$f") (rclone exit $rc)"; failed+=("$(basename "$f")") ;;
    esac
  done
  rclone copy "$MANIFEST" "$DEST/" "${opts[@]}"

  # Verify by hash, not by size: Drive returns MD5, so this is a real integrity
  # check of what landed, and it is the reason `clean` is allowed to delete.
  # Scoped to this group with --include: the staging directory also holds the
  # other group's archives, which are not uploaded yet, and an unscoped check
  # reports those as missing and fails a run that actually succeeded.
  log "verify (md5) $STAGE -> $DEST  [${GROUP}__*.${EXT}]"
  if rclone check "$STAGE" "$DEST" --include "${GROUP}__*.${EXT}" \
       --checkers 2 --log-file "$LOG" --log-level INFO; then
    log "VERIFIED  all staged archives match Drive by MD5"
    touch "$STAGE/.verified_$GROUP"
  else
    log "VERIFY FAILED — do not run clean; re-run upload"
    failed+=("md5-check")
  fi

  if [[ ${#failed[@]} -gt 0 ]]; then
    log "SUMMARY: ${#failed[@]} failure(s): ${failed[*]}"
    return 1
  fi
  log "SUMMARY: group '$GROUP' complete"
}

# --- clean -------------------------------------------------------------------
do_clean () {
  [[ -f "$STAGE/.verified_$GROUP" ]] || {
    echo "group '$GROUP' was not verified against Drive — refusing to delete staged archives" >&2; exit 1; }
  log "clean staged archives for '$GROUP' (verified against Drive)"
  rm -f "$STAGE/${GROUP}__"*.${EXT} "$STAGE/.verified_$GROUP"
  log "clean done; $(df -Pk "$STAGE" | awk 'NR==2{print int($4/1024/1024)}') GiB free"
}

case "$PHASE" in
  tar)    do_tar ;;
  upload) do_upload ;;
  clean)  do_clean ;;
esac
