#!/bin/bash
# fullsim_pp24_full_to_lgd.sh
#
# Take the 6 pp24-fullsim FULL-sample grid tasks from "submitted" to
# "readable by the NTuple processing", WITHOUT downloading ~280 GB to GPFS.
#
#   grid task done -> rucio add-rule to BNL-OSG2_LOCALGROUPDISK -> wait for
#   replication -> symlink farm in ~/usatlasdata/pythia_fullsim_full_sample/
#   -> (optionally) download ONE small slice locally for dev/testing.
#
# Design rationale: Analysis/docs/tracking/pythia_fullsim_pp24_full_sample_skim.md
# §Design Decisions D1.  Measured NTUP size is 30.6 kB/event => the 6 slices are
# ~280 GB; the GPFS data quota cannot hold them.  The future HIJING-overlay full
# sample would be ~20 TB, so read-from-LGD is the architecture regardless.
#
# SAFE TO INTERRUPT AND RE-RUN.  Every stage is idempotent: it checks what is
# already done (rule exists / replica present / symlink correct) and skips it.
#
# Usage (run in tmux; no Claude involved):
#   tmux new -s lgd
#   ./fullsim_pp24_full_to_lgd.sh                 # full chain, polls until done
#   ./fullsim_pp24_full_to_lgd.sh --status        # print state and exit
#   ./fullsim_pp24_full_to_lgd.sh --no-devslice   # skip the local dev-slice download
#
# Requires a valid VOMS proxy with the /atlas/usatlas group:
#   voms-proxy-init -voms atlas -valid 96:00

set -o pipefail

# ---------------------------------------------------------------- configuration
FARM_DIR="/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample"
LOG="${FARM_DIR}/lgd_farm.log"
RULE_FILE="${FARM_DIR}/lgd_rules.txt"          # DID <tab> RULE_ID   (durable state)
RECORD="${FARM_DIR}/merging-record.txt"
LGD_RSE="BNL-OSG2_LOCALGROUPDISK"
SCOPE="user.yuhang"
POLL_MIN=10

# The dev slice kept on GPFS: smallest slice, ~11 GB.  Everything else = LGD only.
DEV_SLICE="pTH70_125"

# jediTaskID  slice   (RE-SKIM 2026-07-21 with the muon_match_L1MU3V branch, VER_TAG=FullJuly2026.v2)
TASKS=(
  "51643327 pTH8_14"
  "51643336 pTH14_24"
  "51643344 pTH24_40"
  "51643353 pTH40_70"
  "51643363 pTH70_125"
  "51643375 pTH125_300"
)

DO_DEVSLICE=1
STATUS_ONLY=0
for a in "$@"; do
  case "$a" in
    --no-devslice) DO_DEVSLICE=0 ;;
    --status)      STATUS_ONLY=1 ;;
    -h|--help)     sed -n '2,/^$/s/^# \{0,1\}//p' "$0"; exit 0 ;;
  esac
done

mkdir -p "$FARM_DIR"
ts()  { date "+%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] $*" | tee -a "$LOG"; }

outds_for() {  # slice -> output dataset name
  echo "${SCOPE}.NTUP.Pythia_5p36TeV_pp_hQCD_DiMu_$1.FullSimPP24.FullJuly2026.v2._EXT0"
}

# ------------------------------------------------------------------ environment
setup_env() {
  source /usatlas/u/yuhanguo/setup.sh >/dev/null 2>&1 || { echo "setup.sh failed"; exit 1; }
  lsetup rucio >/dev/null 2>&1 || { echo "lsetup rucio failed"; exit 1; }
  local left
  left=$(voms-proxy-info --timeleft 2>/dev/null || echo 0)
  if [[ "${left:-0}" -lt 3600 ]]; then
    log "ERROR: VOMS proxy has <1 h left (${left}s). Run: voms-proxy-init -voms atlas -valid 96:00"
    exit 1
  fi
  log "env ready; proxy ${left}s"
}

# --------------------------------------------------------- stage 1: grid status
# echo "done" | "running" | "failed"
grid_state() {
  curl -s --connect-timeout 30 --max-time 60 "https://bigpanda.cern.ch/task/$1/?json" 2>/dev/null |
  python3 -c "
import json,sys
try: d=json.loads(sys.stdin.read())
except: print('running'); raise SystemExit          # transient API error -> keep waiting
t=d.get('task',d); di=t.get('dsinfo',{})
st=t.get('status','unknown')
nf=di.get('nfiles',0) or 0; nff=di.get('nfilesfinished',0) or 0
pct=(nff*100.0/nf) if nf else 0
if st=='done' or (st=='finished' and pct>=90): print('done')
elif st in ('broken','aborted','failed','exhausted'):  print('failed')
else: print('running')
"
}

# --------------------------------------------------- stage 2: replication rule
# Idempotent: reuse an existing rule for (DID, RSE) if one is already there.
ensure_rule() {
  local did="$1" rid
  rid=$(grep -P "^\Q${did}\E\t" "$RULE_FILE" 2>/dev/null | cut -f2 | head -1)
  if [[ -n "$rid" ]]; then echo "$rid"; return 0; fi

  # already has a rule on the RSE from a previous run?
  rid=$(rucio list-rules "${SCOPE}:${did}" 2>/dev/null | awk -v rse="$LGD_RSE" '$0 ~ rse {print $1; exit}')
  if [[ -z "$rid" ]]; then
    rid=$(rucio add-rule "${SCOPE}:${did}" 1 "$LGD_RSE" 2>&1 | grep -oE '^[0-9a-f]{32}' | head -1)
    [[ -z "$rid" ]] && { log "  ERROR: add-rule failed for $did"; return 1; }
    log "  rule created: $rid"
  else
    log "  rule already existed: $rid"
  fi
  printf '%s\t%s\n' "$did" "$rid" >> "$RULE_FILE"
  echo "$rid"
}

rule_ok() {  # rule OK -> 0
  rucio rule-info "$1" 2>/dev/null | grep -qE '^ *State: *OK'
}

# ------------------------------------------------------ stage 3: symlink farm
# Farm names: Pythia_5p36TeV_pp_hQCD_DiMu_<slice>.FullSimPP24.NTUP.partNN.root
# (the NTuple processing TChain-globs "...NTUP.part*.root" -- see D1)
build_farm() {
  local slice="$1" did="$2" n=0 bad=0
  mapfile -t pfns < <(
    rucio list-file-replicas --rse "$LGD_RSE" --pfns "${SCOPE}:${did}" 2>/dev/null |
    grep -E '^(root|gsiftp|davs|/)' | sed -E 's#^[a-z0-9+.-]+://[^/]+##' | sort
  )
  if [[ ${#pfns[@]} -eq 0 ]]; then log "  WARN: no replicas listed for $did"; return 1; fi

  for p in "${pfns[@]}"; do
    n=$((n+1))
    local link
    link=$(printf "%s/Pythia_5p36TeV_pp_hQCD_DiMu_%s.FullSimPP24.NTUP.part%02d.root" "$FARM_DIR" "$slice" "$n")
    if [[ ! -e "$p" ]]; then log "  WARN: pnfs path not visible: $p"; bad=$((bad+1)); continue; fi
    ln -sfn "$p" "$link"
  done
  log "  farm: ${n} symlinks for ${slice} (${bad} unresolved)"
  [[ $bad -eq 0 ]]
}

verify_farm() {  # ROOT must read entries THROUGH the symlinks
  local slice="$1"
  local pat="${FARM_DIR}/Pythia_5p36TeV_pp_hQCD_DiMu_${slice}.FullSimPP24.NTUP.part*.root"
  local n
  n=$(root -l -b -q -e "TChain c(\"HeavyIonD3PD\"); c.Add(\"${pat}\"); printf(\"NENT %lld\n\", c.GetEntries());" 2>/dev/null |
      grep -oP 'NENT \K[0-9]+' | tail -1)
  if [[ -z "$n" || "$n" -eq 0 ]]; then log "  VERIFY FAIL: ${slice} -> 0 entries through farm"; return 1; fi
  log "  VERIFY OK: ${slice} -> ${n} entries through the symlink farm"
  echo "Pythia_5p36TeV_pp_hQCD_DiMu_${slice}.FullSimPP24.NTUP.part*.root | ${LGD_RSE} | $(date +%F) | ${n} entries (LGD symlink farm)" >> "$RECORD"
}

print_status() {
  printf "\n%-12s %-10s %-8s %-6s %s\n" TASK SLICE GRID RULE FARM
  for e in "${TASKS[@]}"; do
    set -- $e; local tid="$1" slice="$2"
    local did; did=$(outds_for "$slice")
    local g r f
    g=$(grid_state "$tid")
    r=$(grep -P "^\Q${did}\E\t" "$RULE_FILE" 2>/dev/null | cut -f2 | head -1)
    if [[ -n "$r" ]]; then rule_ok "$r" && r="OK" || r="repl"; else r="-"; fi
    f=$(ls "${FARM_DIR}/Pythia_5p36TeV_pp_hQCD_DiMu_${slice}.FullSimPP24.NTUP.part"*.root 2>/dev/null | wc -l)
    printf "%-12s %-10s %-8s %-6s %s links\n" "$tid" "$slice" "$g" "$r" "$f"
  done
  echo
}

# ------------------------------------------------------------------------ main
setup_env
touch "$RULE_FILE"

if [[ $STATUS_ONLY -eq 1 ]]; then print_status; exit 0; fi

log "════════ fullsim pp24 FULL sample -> ${LGD_RSE} (farm: ${FARM_DIR}) ════════"

declare -A DONE_FARM
while true; do
  all_done=1
  for e in "${TASKS[@]}"; do
    set -- $e; tid="$1"; slice="$2"
    did=$(outds_for "$slice")

    # already farmed and verified in a previous pass?
    if [[ -n "${DONE_FARM[$slice]:-}" ]]; then continue; fi
    if ls "${FARM_DIR}/Pythia_5p36TeV_pp_hQCD_DiMu_${slice}.FullSimPP24.NTUP.part"*.root >/dev/null 2>&1 \
       && grep -q "DiMu_${slice}\..*symlink farm" "$RECORD" 2>/dev/null; then
      DONE_FARM[$slice]=1; continue
    fi

    g=$(grid_state "$tid")
    case "$g" in
      failed)  log "task $tid ($slice): grid FAILED — needs manual attention"; DONE_FARM[$slice]=1; continue ;;
      running) all_done=0; continue ;;
    esac

    log "task $tid ($slice): grid done -> ensuring LGD rule"
    rid=$(ensure_rule "$did") || { all_done=0; continue; }

    if ! rule_ok "$rid"; then
      log "  rule $rid still replicating (FTS queue can take 1–12 h)"; all_done=0; continue
    fi

    log "  rule OK -> building symlink farm"
    if build_farm "$slice" "$did" && verify_farm "$slice"; then
      DONE_FARM[$slice]=1
      log "  ✅ ${slice} available via LGD farm"
    else
      log "  ${slice} farm incomplete — will retry next pass"; all_done=0
    fi
  done

  [[ $all_done -eq 1 ]] && break
  log "waiting ${POLL_MIN} min…"; sleep $((POLL_MIN*60))
done

# --------------------------------------------- dev slice: one small local copy
# The dev slice lives in its OWN directory, under the CANONICAL name
# "...FullSimPP24.NTUP.root". Two reasons:
#   * the reader only recognises "<base>.NTUP.root" or "<base>.NTUP.part*.root"; a name like
#     ".NTUP.local.root" matches NEITHER and would be silently unreadable;
#   * putting it beside the farm would trip the reader's AMBIGUOUS-input throw (a hadded file
#     AND .part* files for the same slice).
# Use it by pointing the NTuple processing at DEV_DIR via fullsim_input_dir_override.
# NOTE it holds ONE slice, so any run over it needs allow_missing_slices=true and is a
# DIAGNOSTIC run only -- a cross-section-weighted result from a single slice is meaningless.
if [[ $DO_DEVSLICE -eq 1 ]]; then
  DEV_DIR="${FARM_DIR}/local_dev"
  dev_stage="${DEV_DIR}/_download"
  hadded="${DEV_DIR}/Pythia_5p36TeV_pp_hQCD_DiMu_${DEV_SLICE}.FullSimPP24.NTUP.root"
  if [[ -f "$hadded" ]]; then
    log "dev slice already present: $hadded"
  else
    log "downloading dev slice ${DEV_SLICE} (~11 GB) for local iteration"
    mkdir -p "$dev_stage"
    did=$(outds_for "$DEV_SLICE")
    if rucio download --dir "$dev_stage" "${SCOPE}:${did}" >>"$LOG" 2>&1; then
      mapfile -t roots < <(find "$dev_stage" -name '*.root' -type f | sort)
      if [[ ${#roots[@]} -gt 0 ]] && hadd -f "$hadded" "${roots[@]}" >>"$LOG" 2>&1; then
        log "dev slice hadded -> $hadded"
        log "  use: py.fullsim_input_dir_override = \"${DEV_DIR}/\"; py.allow_missing_slices = true;"
        rm -rf "$dev_stage"
      else
        log "ERROR: dev-slice hadd failed (download kept at $dev_stage)"
      fi
    else
      log "ERROR: dev-slice download failed (quota? proxy?)"
    fi
  fi
  # The dev slice needs the SAME AMI dir as the farm (same _pdf production).
  ln -sfn "${FARM_DIR}/ami_info" "${DEV_DIR}/ami_info" 2>/dev/null || true
fi

log "════════ DONE ════════"
print_status
log "NTuple processing: point fullsim_input_dir_override at ${FARM_DIR}"
log "NOTE: requires the TChain-glob change in PythiaAlgCoreT.c InitInputFullsim (D1)."
