#!/bin/bash
# grid_monitor.sh — Monitor grid tasks, download, hadd, validate, update record.
# Run in tmux to survive logout. Renew VOMS proxy from another session on the
# same node when needed (proxy is a shared file at /tmp/x509up_u<uid>).
#
# Usage:
#   ./grid_monitor.sh may2026_skim.txt
#   ./grid_monitor.sh -i 3 may2026_skim.txt          # check every 3 hours
#   ./grid_monitor.sh 50267206 50267236 50267371      # task IDs directly
#
# Input file format: any file containing lines with an 8-digit task ID and
# optionally a user.yuhang..._EXT0 outDS on the same line.  Comment lines
# (#) and decoration are ignored.  Works directly with may2026_skim.txt.

set -o pipefail

# ── Configuration ────────────────────────────────────────────────────────────
INTERVAL_HOURS=5
DATA_BASE="/usatlas/u/yuhanguo/usatlasdata/dimuon_data"
RECORD_FILE="${DATA_BASE}/data-merging-record.txt"
LOG_DIR="${DATA_BASE}"
STATUS_LOG="${LOG_DIR}/grid_monitor_status.log"
ERROR_LOG="${LOG_DIR}/grid_monitor_errors.log"
STATE_FILE="${LOG_DIR}/grid_monitor_state.txt"
BIGPANDA_URL="https://bigpanda.cern.ch/task"
MIN_SUCCESS_PCT=90
# ─────────────────────────────────────────────────────────────────────────────

# ── Parse arguments ──────────────────────────────────────────────────────────
TASK_IDS=()
TASK_OUTDS=()   # parallel array: outDS for each task (may be empty)

while [[ $# -gt 0 ]]; do
	case "$1" in
		-i) INTERVAL_HOURS="$2"; shift 2 ;;
		-h|--help)
			sed -n '2,/^$/s/^# //p' "$0"; exit 0 ;;
		*)
			if [[ -f "$1" ]]; then
				while IFS= read -r line; do
					[[ "$line" =~ ^[[:space:]]*# ]] && continue
					[[ -z "$line" ]] && continue
					tid=$(echo "$line" | grep -oE '\b[0-9]{8}\b' | head -1)
					[[ -z "$tid" ]] && continue
					outds=$(echo "$line" | grep -oE 'user\.yuhang\.[A-Za-z0-9._]+_EXT0' | head -1)
					TASK_IDS+=("$tid")
					TASK_OUTDS+=("${outds:-}")
				done < "$1"
			else
				TASK_IDS+=("$1")
				TASK_OUTDS+=("")
			fi
			shift ;;
	esac
done

if [[ ${#TASK_IDS[@]} -eq 0 ]]; then
	echo "ERROR: no task IDs provided. Usage: $0 [-i HOURS] FILE_OR_TASKIDS..." >&2
	exit 1
fi

INTERVAL=$((INTERVAL_HOURS * 3600))

# ── Logging helpers ──────────────────────────────────────────────────────────
ts() { date "+%Y-%m-%d %H:%M:%S"; }

log_status() {
	local msg="[$(ts)] $*"
	echo "$msg"
	echo "$msg" >> "$STATUS_LOG"
}

log_error() {
	local msg="[$(ts)] ERROR: $*"
	echo "$msg" >&2
	echo "$msg" >> "$ERROR_LOG"
}

# ── State management ─────────────────────────────────────────────────────────
# State file tracks: TASK_ID STATUS (pending|completed|failed)
# On startup, load existing state (supports restart after crash).

declare -A TASK_STATE

load_state() {
	if [[ -f "$STATE_FILE" ]]; then
		while IFS=' ' read -r tid st; do
			[[ -n "$tid" ]] && TASK_STATE["$tid"]="$st"
		done < "$STATE_FILE"
	fi
	for tid in "${TASK_IDS[@]}"; do
		[[ -z "${TASK_STATE[$tid]:-}" ]] && TASK_STATE["$tid"]="pending"
	done
}

save_state() {
	: > "$STATE_FILE"
	for tid in "${TASK_IDS[@]}"; do
		echo "$tid ${TASK_STATE[$tid]}" >> "$STATE_FILE"
	done
}

# ── outDS → directory / filename mapping ─────────────────────────────────────
# user.yuhang.TrigRates.dimuon.PbPb2023data.May2026.v1.part1._EXT0
#   → target_dir=pbpb_2023  output_name=data_pbpb23_part1

map_outds() {
	local outds="$1"
	local dir="" prefix="" part=""

	part=$(echo "$outds" | grep -oP 'part\d+')

	if   [[ "$outds" == *PbPb2023* ]]; then dir="pbpb_2023"; prefix="data_pbpb23"
	elif [[ "$outds" == *PbPb2024* ]]; then dir="pbpb_2024"; prefix="data_pbpb24"
	elif [[ "$outds" == *PbPb2025* ]]; then dir="pbpb_2025"; prefix="data_pbpb25"
	elif [[ "$outds" == *pp2024*   ]]; then dir="pp_2024";   prefix="data_pp24"
	else
		echo "UNKNOWN UNKNOWN"
		return 1
	fi

	echo "${dir} ${prefix}_${part}"
}

# ── Entry counting (from hadd_merge_datasets.sh) ────────────────────────────
count_entries_chain() {
	local tmpList tmpC fn result
	tmpList=$(mktemp /tmp/flist_XXXXXX.txt)
	printf '%s\n' "$@" > "$tmpList"
	tmpC=$(mktemp /tmp/cent_XXXXXX.C)
	fn=$(basename "$tmpC" .C)
	cat > "$tmpC" << 'ROOTEOF'
#include <fstream>
#include <string>
#include <vector>
void FNNAME(const char* listfile) {
    std::ifstream ifs(listfile);
    std::vector<std::string> files;
    std::string line;
    while (std::getline(ifs, line))
        if (!line.empty()) files.push_back(line);
    if (files.empty()) { printf("0\n"); return; }
    TFile* ff = TFile::Open(files[0].c_str());
    std::vector<std::string> trees;
    if (ff && !ff->IsZombie()) {
        TIter nxt(ff->GetListOfKeys()); TKey* k;
        while ((k = (TKey*)nxt()))
            if (strcmp(k->GetClassName(), "TTree") == 0)
                trees.push_back(std::string(k->GetName()));
        ff->Close();
    }
    Long64_t total = 0;
    for (const auto& tn : trees) {
        TChain ch(tn.c_str());
        for (const auto& fp : files) ch.Add(fp.c_str());
        total += ch.GetEntries();
    }
    printf("%lld\n", total);
}
ROOTEOF
	sed -i "s/FNNAME/${fn}/g" "$tmpC"
	result=$(root -l -b -q "${tmpC}(\"${tmpList}\")" 2>/dev/null | grep -E '^[0-9]+$' | tail -1)
	rm -f "$tmpC" "$tmpList"
	echo "${result:-0}"
}

# ── Query BigPanDA API ───────────────────────────────────────────────────────
# Returns: STATUS|SUCCESS_PCT|TASKNAME|ERROR_MSG
query_task() {
	local tid="$1"
	local raw
	raw=$(curl -s --connect-timeout 30 --max-time 60 "${BIGPANDA_URL}/${tid}/?json" 2>/dev/null)
	if [[ -z "$raw" ]]; then
		echo "api_error|0|unknown|curl failed"
		return 1
	fi

	python3 -c "
import json, sys
try:
    d = json.loads(sys.argv[1])
except:
    print('api_error|0|unknown|json parse failed'); sys.exit(0)
task = d.get('task', d)
status = task.get('status', 'unknown')
taskname = task.get('taskname', '')
errmsg = (task.get('errordialog', '') or '')[:200]
di = task.get('dsinfo', {})
nf = di.get('nfiles', 0) or 0
nff = di.get('nfilesfinished', 0) or 0
pct = (nff * 100.0 / nf) if nf > 0 else 0.0
print(f'{status}|{pct:.1f}|{taskname}|{errmsg}')
" "$raw"
}

# ── Check VOMS proxy ────────────────────────────────────────────────────────
check_proxy() {
	local timeleft
	timeleft=$(voms-proxy-info --timeleft 2>/dev/null || echo 0)
	if [[ "$timeleft" -lt 600 ]]; then
		log_error "VOMS proxy has <10 min remaining (${timeleft}s). Renew from another session on this node: voms-proxy-init -voms atlas"
		return 1
	fi
	return 0
}

# ── Download + hadd + validate for one task ──────────────────────────────────
process_task() {
	local tid="$1"
	local outds="$2"

	local mapping target_subdir output_name
	mapping=$(map_outds "$outds")
	if [[ "$mapping" == "UNKNOWN UNKNOWN" ]]; then
		log_error "task $tid: cannot map outDS '$outds' to directory/filename"
		return 1
	fi
	target_subdir=$(echo "$mapping" | cut -d' ' -f1)
	output_name=$(echo "$mapping" | cut -d' ' -f2)

	local target_dir="${DATA_BASE}/${target_subdir}"
	local output_file="${target_dir}/${output_name}.root"
	local rucio_did="user.yuhang:${outds}"
	local rucio_subdir="${target_dir}/${outds}"

	mkdir -p "$target_dir"

	# Back up existing output file if present
	if [[ -f "$output_file" ]]; then
		local bak="${output_file%.root}.bak_$(date +%Y%m%d).root"
		log_status "task $tid: backing up existing $output_name.root → $(basename $bak)"
		mv "$output_file" "$bak"
	fi

	# Download
	log_status "task $tid: downloading $outds → $target_subdir/"
	if ! rucio download --dir "$target_dir" "$rucio_did" 2>&1 | tee -a "$STATUS_LOG"; then
		log_error "task $tid: rucio download failed for $outds"
		return 1
	fi

	# Verify download directory exists and has ROOT files
	if [[ ! -d "$rucio_subdir" ]]; then
		log_error "task $tid: expected download dir $rucio_subdir not found"
		return 1
	fi
	mapfile -t roots < <(find "$rucio_subdir" -maxdepth 1 -name "*.root" -type f | sort)
	if [[ ${#roots[@]} -eq 0 ]]; then
		log_error "task $tid: no ROOT files in $rucio_subdir"
		return 1
	fi
	log_status "task $tid: downloaded ${#roots[@]} files"

	# Count entries before merge
	local pre_count
	pre_count=$(count_entries_chain "${roots[@]}")
	log_status "task $tid: pre-merge entries = $pre_count"

	# hadd
	log_status "task $tid: hadd (${#roots[@]} files) → $output_name.root"
	if ! hadd -f "$output_file" "${roots[@]}" >> "$STATUS_LOG" 2>&1; then
		log_error "task $tid: hadd failed for $output_name"
		rm -f "$output_file"
		return 1
	fi

	# Count entries after merge
	local post_count
	post_count=$(count_entries_chain "$output_file")
	log_status "task $tid: post-merge entries = $post_count"

	# Validate
	if [[ "$post_count" -ne "$pre_count" ]]; then
		log_error "task $tid: ENTRY MISMATCH for $output_name: pre=$pre_count post=$post_count"
		echo "${output_name}.root | ${outds} | $(date +%F) | ${post_count} entries !!! MISMATCH: source=${pre_count} merged=${post_count}" >> "$RECORD_FILE"
		return 1
	fi

	# Record and cleanup
	echo "${output_name}.root | ${outds} | $(date +%F) | ${post_count} entries" >> "$RECORD_FILE"
	log_status "task $tid: validated OK — ${post_count} entries. Cleaning up download dir."
	rm -rf "$rucio_subdir"
	return 0
}

# ── Resolve outDS for a task (query API if not provided) ─────────────────────
resolve_outds() {
	local tid="$1"
	local known_outds="$2"

	if [[ -n "$known_outds" ]]; then
		echo "$known_outds"
		return 0
	fi

	local raw
	raw=$(curl -s --connect-timeout 30 --max-time 60 "${BIGPANDA_URL}/${tid}/?json" 2>/dev/null)
	python3 -c "
import json, sys
try:
    d = json.loads(sys.argv[1])
except:
    sys.exit(1)
datasets = d.get('datasets', [])
for dd in datasets:
    if dd.get('type') == 'output':
        name = dd.get('datasetname', '')
        # Strip trailing JEDI suffixes (e.g. .663585806.663585806)
        # Keep up to _EXT0
        idx = name.find('_EXT0')
        if idx >= 0:
            print(name[:idx+5])
            sys.exit(0)
# Fallback: derive from taskname
task = d.get('task', d)
tn = task.get('taskname', '').rstrip('/')
if tn:
    print(tn + '_EXT0')
" "$raw"
}

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

# Environment setup (ROOT for hadd, rucio for download)
log_status "Setting up environment..."
source ~/setup.sh || { echo "ERROR: setup.sh failed" >&2; exit 1; }
lsetup rucio      || { echo "ERROR: lsetup rucio failed" >&2; exit 1; }
log_status "Environment ready. root=$(command -v root) rucio=$(command -v rucio)"

load_state

log_status "════════════════════════════════════════════════════════════════"
log_status "Grid monitor started: ${#TASK_IDS[@]} tasks, checking every ${INTERVAL_HOURS}h"
log_status "  Status log: $STATUS_LOG"
log_status "  Error log:  $ERROR_LOG"
log_status "  State file: $STATE_FILE"
log_status "════════════════════════════════════════════════════════════════"

for i in "${!TASK_IDS[@]}"; do
	log_status "  ${TASK_IDS[$i]}  ${TASK_OUTDS[$i]:-<will resolve>}  [${TASK_STATE[${TASK_IDS[$i]}]}]"
done

cycle=0
while true; do
	((cycle++))
	log_status ""
	log_status "━━━━━━━━━━ Cycle $cycle  $(ts) ━━━━━━━━━━"

	pending=0
	for i in "${!TASK_IDS[@]}"; do
		tid="${TASK_IDS[$i]}"
		[[ "${TASK_STATE[$tid]}" != "pending" ]] && continue
		((pending++))

		log_status "Checking task $tid ..."
		result=$(query_task "$tid")
		IFS='|' read -r status pct taskname errmsg <<< "$result"

		log_status "  status=$status  success=${pct}%  taskname=$taskname"
		[[ -n "$errmsg" ]] && log_status "  errordialog: $errmsg"

		case "$status" in
			done)
				log_status "  Task $tid DONE (100%). Proceeding to download+merge."
				;;
			finished)
				if (( $(echo "$pct >= $MIN_SUCCESS_PCT" | bc -l) )); then
					log_status "  Task $tid FINISHED (${pct}% >= ${MIN_SUCCESS_PCT}%). Proceeding to download+merge."
				else
					log_error "task $tid: FINISHED with low success rate (${pct}% < ${MIN_SUCCESS_PCT}%). Marking failed."
					TASK_STATE["$tid"]="failed"
					save_state
					continue
				fi
				;;
			broken|aborted|failed|exhausted)
				log_error "task $tid: terminal status '$status'. $errmsg"
				TASK_STATE["$tid"]="failed"
				save_state
				continue
				;;
			api_error)
				log_error "task $tid: API query failed ($errmsg). Will retry next cycle."
				continue
				;;
			*)
				log_status "  Task $tid still in progress ($status). Skipping."
				continue
				;;
		esac

		# Task is ready for download — check proxy first
		if ! check_proxy; then
			log_status "  Skipping download (proxy expired). Will retry next cycle."
			continue
		fi

		# Resolve outDS
		outds=$(resolve_outds "$tid" "${TASK_OUTDS[$i]:-}")
		if [[ -z "$outds" ]]; then
			log_error "task $tid: could not resolve output dataset name"
			continue
		fi
		TASK_OUTDS[$i]="$outds"
		log_status "  outDS: $outds"

		# Download, hadd, validate
		if process_task "$tid" "$outds"; then
			TASK_STATE["$tid"]="completed"
			log_status "  Task $tid: COMPLETED successfully."
		else
			log_error "task $tid: download/merge failed. Will NOT retry automatically — investigate."
			TASK_STATE["$tid"]="failed"
		fi
		save_state
	done

	# Recount pending
	pending=0
	for tid in "${TASK_IDS[@]}"; do
		[[ "${TASK_STATE[$tid]}" == "pending" ]] && ((pending++))
	done

	# Summary
	completed=0; failed=0
	for tid in "${TASK_IDS[@]}"; do
		[[ "${TASK_STATE[$tid]}" == "completed" ]] && ((completed++))
		[[ "${TASK_STATE[$tid]}" == "failed" ]] && ((failed++))
	done
	log_status "Cycle $cycle done: $completed completed, $failed failed, $pending pending"

	if [[ $pending -eq 0 ]]; then
		log_status ""
		log_status "════════════════════════════════════════════════════════════════"
		log_status "All tasks resolved. $completed completed, $failed failed."
		log_status "════════════════════════════════════════════════════════════════"
		exit 0
	fi

	log_status "Sleeping ${INTERVAL_HOURS}h until next check ($(date -d "+${INTERVAL} seconds" "+%Y-%m-%d %H:%M"))..."
	sleep "$INTERVAL"
done
