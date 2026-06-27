#!/usr/bin/env bash
# Copy/sync iacobus CHIME_DSA_Codetections → Google Drive (jakobtfaber@gmail.com).
# Runs rclone ON iacobus (direct upload; jakob-mbp is orchestrator only).
#
# Prerequisite — gdrive remote (once per host that runs rclone):
#   rclone config   # name: gdrive-jakob, type: drive, scope: drive, account: jakobtfaber@gmail.com
# OAuth needs a browser. Headless iacobus flow:
#   1. On jakob-mbp:  rclone authorize "drive"
#   2. Copy the {"access_token":...} blob
#   3. On iacobus:    rclone config  → n) New remote → gdrive-jakob → paste token at "config_token"
# Verify: rclone about gdrive-jakob: && rclone lsd gdrive-jakob:Research/
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
MANIFEST="${MANIFEST:-$REPO_ROOT/codetections_manifest.yaml}"
IACOBUS="${IACOBUS_HOST:-iacobus}"
REMOTE="${RCLONE_REMOTE:-gdrive-jakob}"
GDRIVE_ROOT="${GDRIVE_ROOT:-Research/CHIME_DSA_Codetections}"
SRC="${IACOBUS_SRC:-/Users/iacobus/Research/CHIME_DSA_Codetections}"
LOG_DIR="${LOG_DIR:-$HOME/logs/gdrive_transfers}"
RCLONE_OPTS=(--progress --stats-one-line --stats 30s -v)
TRANSFER="${TRANSFER:-copy}" # copy | sync

usage() {
  cat <<EOF >&2
Usage: $0 [--dry-run] [--subdir NAME] [--verify-only] [--sync]

  --dry-run       rclone dry-run (no bytes moved)
  --subdir NAME   transfer one top-level subdir only (e.g. metadata)
  --verify-only   sentinel SHA-256 check on remote; no transfer
  --sync          use rclone sync instead of copy (destructive on remote extras)

Env: RCLONE_REMOTE (default gdrive-jakob), GDRIVE_ROOT, IACOBUS_HOST, LOG_DIR
EOF
  exit 1
}

DRY=0
SUBDIR=""
VERIFY_ONLY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY=1; shift ;;
    --subdir) SUBDIR="$2"; shift 2 ;;
    --verify-only) VERIFY_ONLY=1; shift ;;
    --sync) TRANSFER=sync; shift ;;
    -h|--help) usage ;;
    *) usage ;;
  esac
done

NAME="iacobus_to_gdrive"
[[ -n "$SUBDIR" ]] && NAME="${NAME}_${SUBDIR}"
LOG="$LOG_DIR/${NAME}.log"
mkdir -p "$LOG_DIR"

LOCAL_SRC="$SRC/"
REMOTE_DST="${REMOTE}:${GDRIVE_ROOT}/"
[[ -n "$SUBDIR" ]] && LOCAL_SRC="${SRC}/${SUBDIR}/" && REMOTE_DST="${REMOTE}:${GDRIVE_ROOT}/${SUBDIR}/"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

sentinel_verify() {
  python3 - "$MANIFEST" "$REMOTE" "$GDRIVE_ROOT" <<'PY' | ssh -o BatchMode=yes "$IACOBUS" "bash -s" >>"$LOG" 2>&1
import hashlib, json, subprocess, sys, yaml

manifest_path, remote, gdrive_root = sys.argv[1:4]
manifest = yaml.safe_load(open(manifest_path))
meta = manifest["subdirs"]["metadata"]
rel = meta["sentinel_path"]
expected = meta["sentinel_sha256"]
nbytes = meta["sentinel_sha256_prefix_bytes"]
remote_path = f"{remote}:{gdrive_root}/{rel}"

proc = subprocess.run(
    ["rclone", "cat", remote_path],
    capture_output=True,
)
if proc.returncode != 0:
    print(f"SENTINEL FAIL: rclone cat {remote_path}: {proc.stderr.decode().strip()}")
    sys.exit(1)
data = proc.stdout[:nbytes]
got = hashlib.sha256(data).hexdigest()
if got != expected:
    print(f"SENTINEL FAIL: {rel} sha256 prefix mismatch")
    print(f"  expected: {expected}")
    print(f"  got:      {got}")
    sys.exit(1)
print(f"SENTINEL PASS: {rel} ({nbytes} B prefix sha256 match)")
PY
}

if (( VERIFY_ONLY )); then
  log "VERIFY-ONLY sentinel check remote=$REMOTE dst=$GDRIVE_ROOT"
  sentinel_verify
  echo "OK sentinel (log: $LOG)"
  exit 0
fi

RCLONE_DRY=()
(( DRY )) && RCLONE_DRY=(--dry-run)

log "START $NAME transfer=$TRANSFER dry_run=$DRY"
log "  src: iacobus:$LOCAL_SRC"
log "  dst: $REMOTE_DST"

# Preflight: remote exists on iacobus
if ! ssh -o BatchMode=yes "$IACOBUS" "rclone listremotes" 2>>"$LOG" | grep -q "^${REMOTE}:$"; then
  log "FAILED preflight: remote ${REMOTE}: not configured on $IACOBUS"
  echo "FAIL — configure gdrive remote on iacobus first (see script header). Log: $LOG" >&2
  exit 1
fi

RSYNC_CMD=(rclone "$TRANSFER" "${RCLONE_DRY[@]}" "${RCLONE_OPTS[@]}" "$LOCAL_SRC" "$REMOTE_DST")
if ssh -o BatchMode=yes "$IACOBUS" "${RSYNC_CMD[@]}" >>"$LOG" 2>&1; then
  log "SUCCESS rclone $TRANSFER"
else
  log "FAILED rclone $TRANSFER"
  echo "FAIL $NAME — see $LOG" >&2
  exit 1
fi

if (( ! DRY )) && [[ -z "$SUBDIR" || "$SUBDIR" == "metadata" ]]; then
  log "sentinel verify metadata/interveners.ipynb"
  sentinel_verify || { echo "FAIL sentinel — see $LOG" >&2; exit 1; }
fi

echo "OK $NAME (log: $LOG)"
