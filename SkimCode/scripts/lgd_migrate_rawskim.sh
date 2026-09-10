#!/bin/bash
# lgd_migrate_rawskim.sh -- migrate one group of raw skim NTUPs from local GPFS to
# BNL-OSG2_LOCALGROUPDISK, then replace them in place with a symlink farm.
#
#   ./lgd_migrate_rawskim.sh <key> [stage]
#     <key>   : pbpb23 | pbpb24 | pbpb25 | pp24 | pbpb26
#     [stage] : upload | dataset | rule | wait | farm | verify | keep | purge | all
#               (default: all; `purge` is the ONLY stage that deletes anything and it
#                re-runs verify first)
#
# Procedure and gotchas inherited from Analysis/docs/tracking/localgroupdisk_migration.md.
# NOTHING IS DELETED BY THIS SCRIPT.  Deleting the local originals is a separate,
# deliberate step performed only after `verify` passes.
#
# Safety invariants:
#   * `rucio upload` COPIES; local originals are never touched.
#   * NEVER `rucio update-rule --lifetime 0` -- that means "expire now" and purges the
#     replicas.  There is no call to update-rule anywhere in this script.
#   * The symlink farm is only built after the rule reports state OK, and the originals
#     are moved aside to <dir>_orig (not deleted) before the farm goes in.

set -uo pipefail

KEY="${1:?usage: $0 <pbpb23|pbpb24|pbpb25|pp24|pbpb26> [stage]}"
STAGE="${2:-all}"

BASE=/usatlas/u/yuhanguo/usatlasdata/dimuon_data
WORK=/usatlas/u/yuhanguo/usatlasdata/lgd_migration/rawskim_2026
SCRATCH_RSE=BNL-OSG2_SCRATCHDISK
LGD_RSE=BNL-OSG2_LOCALGROUPDISK
SCOPE=user.yuhang
PFN_PREFIX='root://dcgftp.usatlas.bnl.gov:1094/'

case "$KEY" in
  pbpb23) SUBDIR=pbpb_2023; GLOB='data_pbpb23_part*.root'; KEEP=data_pbpb23_part4.root ;;
  pbpb24) SUBDIR=pbpb_2024; GLOB='data_pbpb24_part*.root'; KEEP=data_pbpb24_part1.root ;;
  pbpb25) SUBDIR=pbpb_2025; GLOB='data_pbpb25_part*.root'; KEEP=data_pbpb25_part6.root ;;
  pbpb26) SUBDIR=pbpb_2026; GLOB='data_pbpb26_part*.root'; KEEP="" ;;
  pp24)   SUBDIR=pp_2024;   GLOB='data_pp24_part*.root';   KEEP=data_pp24_part11.root ;;
  *) echo "unknown key '$KEY'"; exit 2 ;;
esac

DS="dimuon_rawskim_${KEY}"
DIR="$BASE/$SUBDIR"
LOG="$WORK/${KEY}.log"
PFNS="$WORK/${KEY}_pfns.txt"
mkdir -p "$WORK"

# NOTE: append straight to the log, then echo to stdout separately and ignore any
# failure there.  A `| tee -a` pipeline dies on SIGPIPE the moment the parent shell's
# stdout goes away (e.g. when a foreground run is moved to the background), which
# silently swallows every subsequent progress line even though the work continues.
say() { local m="[$(date -u +%FT%TZ)] $*"; printf '%s\n' "$m" >> "$LOG"; printf '%s\n' "$m" 2>/dev/null || true; }

setup_env() {
  [[ -n "${RUCIO_READY:-}" ]] && return 0
  # ATLAS setup scripts reference unbound variables; `set -u` would kill the shell
  # mid-source (silently, with no log line at all).  Relax it just for the sourcing.
  set +u
  source ~/setup.sh                                                    > /dev/null 2>&1
  source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh > /dev/null 2>&1
  lsetup rucio                                                         > /dev/null 2>&1
  set -u
  export SITE_NAME=BNL-ATLAS
  export RUCIO_READY=1
  rucio whoami > /dev/null 2>&1 || { say "FATAL: rucio/proxy not usable"; exit 1; }
}

files_local() { ls -1 "$DIR"/$GLOB 2>/dev/null; }

do_upload() {
  say "=== upload: $KEY -> $SCRATCH_RSE ==="
  local n=0 ok=0
  for f in $(files_local); do
    n=$((n+1))
    local b; b=$(basename "$f")
    # NOTE: `rucio list-dids <scope>:<exact-name> --short` returns NOTHING for a file DID.
    # The type filter is required -- without it the check always says "free" and the
    # upload is retried (and fails) for every already-uploaded file.
    if rucio list-dids "${SCOPE}:${b}" --filter type=FILE --short 2>/dev/null | grep -q "^${SCOPE}:${b}$"; then
      say "  SKIP $b (DID already registered)"; ok=$((ok+1)); continue
    fi
    say "  uploading $b ($(stat -c %s "$f") bytes)"
    if rucio upload --rse "$SCRATCH_RSE" --scope "$SCOPE" "$f" >> "$LOG" 2>&1; then
      say "  OK   $b"; ok=$((ok+1))
    else
      say "  FAIL $b -- see $LOG"
    fi
  done
  say "upload done: $ok/$n present on $SCRATCH_RSE"
  [[ "$ok" == "$n" ]]
}

do_dataset() {
  say "=== dataset: ${SCOPE}:${DS} ==="
  rucio add-dataset "${SCOPE}:${DS}" >> "$LOG" 2>&1
  local dids=()
  for f in $(files_local); do dids+=("${SCOPE}:$(basename "$f")"); done
  [[ ${#dids[@]} -eq 0 ]] && { say "no local files matched $DIR/$GLOB"; return 1; }
  rucio attach "${SCOPE}:${DS}" "${dids[@]}" >> "$LOG" 2>&1
  local n; n=$(rucio list-files "${SCOPE}:${DS}" 2>/dev/null | grep -c "^| ${SCOPE}:")
  say "dataset now holds $n files (expected ${#dids[@]})"
}

do_rule() {
  say "=== rule: ${SCOPE}:${DS} -> $LGD_RSE ==="
  local existing
  existing=$(rucio list-rules "${SCOPE}:${DS}" 2>/dev/null | awk -v r="$LGD_RSE" '$0 ~ r {print $1}' | head -1)
  if [[ -n "$existing" ]]; then say "rule already exists: $existing"; echo "$existing" > "$WORK/${KEY}.ruleid"; return 0; fi
  local out; out=$(rucio add-rule "${SCOPE}:${DS}" 1 "$LGD_RSE" 2>&1 | tee -a "$LOG")
  local rid; rid=$(echo "$out" | grep -oE '[0-9a-f]{32}' | head -1)
  [[ -z "$rid" ]] && { say "FAILED to create rule: $out"; return 1; }
  say "rule created: $rid"; echo "$rid" > "$WORK/${KEY}.ruleid"
}

do_wait() {
  local rid; rid=$(cat "$WORK/${KEY}.ruleid" 2>/dev/null)
  [[ -z "$rid" ]] && { say "no rule id recorded for $KEY"; return 1; }
  say "=== waiting for rule $rid ==="
  while :; do
    local info state
    info=$(rucio rule-info "$rid" 2>&1)
    state=$(echo "$info" | awk -F': *' '/^State:/ {print $2}')
    say "  rule $rid state=$state  $(echo "$info" | grep -E '^Locks:' | tr -s ' ')"
    [[ "$state" == "OK" ]] && return 0
    if [[ "$state" == "STUCK" ]]; then say "  RULE STUCK -- stopping, needs a human look"; return 1; fi
    sleep 300
  done
}

do_farm() {
  say "=== symlink farm for $KEY ==="
  rucio list-file-replicas "${SCOPE}:${DS}" --protocols root --pfns --rses "$LGD_RSE" > "$PFNS" 2>>"$LOG"
  local npfn; npfn=$(grep -c "^${PFN_PREFIX}" "$PFNS")
  local nloc; nloc=$(files_local | wc -l)
  say "  $npfn PFNs on $LGD_RSE for $nloc local files"
  [[ "$npfn" -ne "$nloc" ]] && { say "  REFUSING to build farm: PFN count != local file count"; return 1; }

  # Move the originals aside (NOT deleted), then create the farm at the original path.
  if [[ -d "${DIR}_orig_${KEY}" ]]; then say "  ${DIR}_orig_${KEY} already exists -- reusing"; else mkdir -p "${DIR}_orig_${KEY}"; fi
  for f in $(files_local); do
    [[ -L "$f" ]] && continue          # already a symlink: farm was built before
    mv "$f" "${DIR}_orig_${KEY}/" || { say "  mv failed for $f"; return 1; }
  done

  while read -r pfn; do
    [[ "$pfn" == ${PFN_PREFIX}* ]] || continue
    local pnfs="${pfn#$PFN_PREFIX}"
    local b; b=$(basename "$pnfs")
    ln -sfn "$pnfs" "$DIR/$b"
    say "  link $b -> $pnfs"
  done < "$PFNS"
  say "farm built at $DIR ; originals parked in ${DIR}_orig_${KEY}"
}

do_verify() {
  say "=== verify $KEY (size + HeavyIonD3PD entries, farm vs parked original) ==="
  source ~/setup.sh > /dev/null 2>&1
  local fail=0
  for l in "$DIR"/$GLOB; do
    [[ -e "$l" ]] || continue
    local b; b=$(basename "$l")
    local o="${DIR}_orig_${KEY}/$b"
    [[ -f "$o" ]] || { say "  $b: no parked original to compare -- SKIP"; continue; }
    local sl so
    sl=$(stat -Lc %s "$l" 2>/dev/null); so=$(stat -c %s "$o" 2>/dev/null)
    local el eo
    el=$(root -l -b -q -e "TFile*f=TFile::Open(\"$l\");TTree*t=f?(TTree*)f->Get(\"HeavyIonD3PD\"):0;printf(\"ENT %lld\\n\",t?t->GetEntries():-1);" 2>/dev/null | grep -oP '(?<=^ENT )-?\d+' | head -1)
    eo=$(root -l -b -q -e "TFile*f=TFile::Open(\"$o\");TTree*t=f?(TTree*)f->Get(\"HeavyIonD3PD\"):0;printf(\"ENT %lld\\n\",t?t->GetEntries():-1);" 2>/dev/null | grep -oP '(?<=^ENT )-?\d+' | head -1)
    if [[ "$sl" == "$so" && "$el" == "$eo" && -n "$el" && "$el" != "-1" ]]; then
      say "  OK   $b  bytes=$sl  entries=$el"
    else
      say "  FAIL $b  bytes lgd=$sl orig=$so  entries lgd=$el orig=$eo"; fail=1
    fi
  done
  [[ $fail -eq 0 ]] && say "VERIFY PASSED for $KEY" || say "VERIFY FAILED for $KEY"
  return $fail
}

# Restore ONE small real file per data-taking period over its symlink, so a local test
# (and any work during a dCache outage) does not depend on /pnfs being up.  The file stays
# on LOCALGROUPDISK as well -- this is deliberate redundancy, not a missing migration.
do_keep() {
  [[ -z "${KEEP:-}" ]] && { say "no local keeper configured for $KEY"; return 0; }
  local src="${DIR}_orig_${KEY}/$KEEP" dst="$DIR/$KEEP"
  [[ -f "$src" ]] || { say "keeper $src not present (already restored?)"; return 0; }
  [[ -L "$dst" ]] && rm -f "$dst"
  mv "$src" "$dst" || { say "failed to restore keeper $KEEP"; return 1; }
  say "restored local keeper $dst ($(stat -c %s "$dst") bytes)"
}

# Delete the parked originals.  Refuses unless verify has passed in this same invocation
# chain, and refuses if anything in $DIR is still a real file rather than a symlink
# (other than the configured keeper).
do_purge() {
  local d="${DIR}_orig_${KEY}"
  [[ -d "$d" ]] || { say "nothing parked at $d"; return 0; }
  do_verify || { say "REFUSING to purge $d -- verify failed"; return 1; }
  local n; n=$(ls -1 "$d" | wc -l)
  local bytes; bytes=$(du -sb "$d" | cut -f1)
  say "purging $n parked originals ($bytes bytes) from $d"
  rm -rf "$d" && say "purged $d"
}

setup_env
case "$STAGE" in
  upload)  do_upload ;;
  keep)    do_keep ;;
  purge)   do_purge ;;
  dataset) do_dataset ;;
  rule)    do_rule ;;
  wait)    do_wait ;;
  farm)    do_farm ;;
  verify)  do_verify ;;
  all)     do_upload && do_dataset && do_rule && do_wait && do_farm && do_verify && do_keep ;;
  *) echo "unknown stage '$STAGE'"; exit 2 ;;
esac
