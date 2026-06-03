#!/bin/bash
# grid_monitor.sh — Monitor grid tasks, download, hadd, validate, update record.
# Supports parallel execution across multiple nodes via flock on the shared
# state file (GPFS). Run one instance per node in tmux.
#
# Usage:
#   ./grid_monitor.sh may2026_skim.txt                # single node, default 10min poll
#   ./grid_monitor.sh -i 15 may2026_skim.txt          # poll every 15 min when idle
#   ./grid_monitor.sh 50267206 50267236 50267371       # task IDs directly
#
# Multi-node: run the same command on each node (e.g. one tmux pane per node).
# Each instance claims one task at a time; flock prevents double-claiming.
# Renew VOMS proxy from another session when needed.
#
# Input file format: any file containing lines with an 8-digit task ID and
# optionally a user.yuhang..._EXT0 outDS on the same line.  Comment lines
# (#) and decoration are ignored.  Works directly with may2026_skim.txt.
#
# State file format: TASK_ID STATE [HOSTNAME EPOCH_TIMESTAMP]
#   States: pending → ready → downloading → completed | failed

set -o pipefail

# ── Configuration ────────────────────────────────────────────────────────────
POLL_INTERVAL_MIN=10
STALE_TIMEOUT=10800   # 3 hours: reclaim downloads from crashed workers
HADD_CHUNK_MAX_FILES=100  # max files per chunk in chunked fallback
ANALYSIS_CODE_DIR="/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/NTupleProcessingCode"
DATA_BASE="/usatlas/u/yuhanguo/usatlasdata/dimuon_data"
RECORD_FILE="${DATA_BASE}/data-merging-record.txt"
LOG_DIR="${DATA_BASE}"
STATUS_LOG="${LOG_DIR}/grid_monitor_status.log"
ERROR_LOG="${LOG_DIR}/grid_monitor_errors.log"
STATE_FILE="${LOG_DIR}/grid_monitor_state.txt"
LOCK_FILE="${STATE_FILE}.lock"
BIGPANDA_URL="https://bigpanda.cern.ch/task"
MIN_SUCCESS_PCT=90
MY_HOSTNAME=$(hostname -s)
# ─────────────────────────────────────────────────────────────────────────────

# ── Parse arguments ──────────────────────────────────────────────────────────
TASK_IDS=()
TASK_OUTDS=()   # parallel array: outDS for each task (may be empty)

while [[ $# -gt 0 ]]; do
	case "$1" in
		-i) POLL_INTERVAL_MIN="$2"; shift 2 ;;
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
	echo "ERROR: no task IDs provided. Usage: $0 [-i POLL_MINUTES] FILE_OR_TASKIDS..." >&2
	exit 1
fi

POLL_INTERVAL=$((POLL_INTERVAL_MIN * 60))

# Build associative array for outDS lookup by task ID
declare -A OUTDS_MAP
for i in "${!TASK_IDS[@]}"; do
	[[ -n "${TASK_OUTDS[$i]}" ]] && OUTDS_MAP["${TASK_IDS[$i]}"]="${TASK_OUTDS[$i]}"
done

# ── Logging helpers ──────────────────────────────────────────────────────────
ts() { date "+%Y-%m-%d %H:%M:%S"; }

log_status() {
	local msg="[$(ts)] [$MY_HOSTNAME] $*"
	echo "$msg"
	echo "$msg" >> "$STATUS_LOG"
}

log_error() {
	local msg="[$(ts)] [$MY_HOSTNAME] ERROR: $*"
	echo "$msg" >&2
	echo "$msg" >> "$ERROR_LOG"
}

# ── Locking helpers ──────────────────────────────────────────────────────────
# fd 200 is opened once; flock/unlock operate on it across all functions.
exec 200>"$LOCK_FILE"
lock_state() { flock -w 30 200; }
unlock_state() { flock -u 200; }

# ── State management ─────────────────────────────────────────────────────────

init_state() {
	lock_state
	[[ -f "$STATE_FILE" ]] || touch "$STATE_FILE"
	for tid in "${TASK_IDS[@]}"; do
		if ! grep -q "^${tid} " "$STATE_FILE"; then
			echo "$tid pending" >> "$STATE_FILE"
		fi
	done
	unlock_state
}

# Compare-and-swap: only update if current state starts with expected_prefix
set_task_state_if() {
	local tid="$1" expected_prefix="$2" new_state="$3"
	lock_state
	if grep -q "^${tid} ${expected_prefix}" "$STATE_FILE"; then
		sed -i "s/^${tid} ${expected_prefix}.*/${tid} ${new_state}/" "$STATE_FILE"
	fi
	unlock_state
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

	echo "$raw" | python3 -c "
import json, sys
try:
    d = json.loads(sys.stdin.read())
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
"
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

# ── Recursive chunked hadd ──────────────────────────────────────────────────
# Merges an array of ROOT files into one output. If direct hadd fails,
# recursively splits into halves until chunks are small enough (≤50 files).
RECURSIVE_HADD_CHUNK=50
recursive_hadd() {
	local output="$1"; shift
	local files=("$@")
	local n=${#files[@]}
	local label="${output##*/}"

	# Small enough to merge directly (no fallback)
	if (( n <= RECURSIVE_HADD_CHUNK )); then
		hadd -f "$output" "${files[@]}" >> "$STATUS_LOG" 2>&1
		return $?
	fi

	# Try direct hadd first
	if hadd -f "$output" "${files[@]}" >> "$STATUS_LOG" 2>&1; then
		return 0
	fi

	# Direct hadd failed — split and recurse
	rm -f "$output"
	log_status "  recursive_hadd: splitting $n files for $label"
	local mid=$(( n / 2 ))
	local tmp_a="${output%.root}_chunka.root"
	local tmp_b="${output%.root}_chunkb.root"

	if ! recursive_hadd "$tmp_a" "${files[@]:0:$mid}"; then
		rm -f "$tmp_a"
		return 1
	fi
	if ! recursive_hadd "$tmp_b" "${files[@]:$mid}"; then
		rm -f "$tmp_a" "$tmp_b"
		return 1
	fi

	if ! hadd -f "$output" "$tmp_a" "$tmp_b" >> "$STATUS_LOG" 2>&1; then
		log_error "  recursive_hadd: final merge failed for $label"
		rm -f "$tmp_a" "$tmp_b" "$output"
		return 1
	fi
	rm -f "$tmp_a" "$tmp_b"
	return 0
}

# ── Code-to-dataset mapping (for auto-update after chunking) ────────────────
# Returns: EXTRAS_FILE SUB_FILE RUN_YEAR_KEY
get_code_update_info() {
	local target_subdir="$1"
	case "$target_subdir" in
		pp_2024)   echo "$ANALYSIS_CODE_DIR/PPExtras.c   $ANALYSIS_CODE_DIR/run_pp_24.sub   24" ;;
		pbpb_2023) echo "$ANALYSIS_CODE_DIR/PbPbExtras.c $ANALYSIS_CODE_DIR/run_pbpb_23.sub 23" ;;
		pbpb_2024) echo "$ANALYSIS_CODE_DIR/PbPbExtras.c $ANALYSIS_CODE_DIR/run_pbpb_24.sub 24" ;;
		pbpb_2025) echo "$ANALYSIS_CODE_DIR/PbPbExtras.c $ANALYSIS_CODE_DIR/run_pbpb_25.sub 25" ;;
		*) echo "" ;;
	esac
}

# Update file_batch_max in C++ source and queue count in .sub, then git commit.
update_source_for_new_max() {
	local extras_file="$1" sub_file="$2" run_year_key="$3" new_max="$4"

	if [[ ! -f "$extras_file" ]]; then
		log_error "Cannot update source: $extras_file not found"
		return 1
	fi

	sed -i "s/{${run_year_key}, [0-9]*}/{${run_year_key}, ${new_max}}/" "$extras_file"
	if grep -q "{${run_year_key}, ${new_max}}" "$extras_file"; then
		log_status "  Updated $extras_file: {$run_year_key, $new_max} — verified"
	else
		log_error "  FAILED to update $extras_file: {$run_year_key, $new_max} not found after sed"
		return 1
	fi

	if [[ -f "$sub_file" ]]; then
		sed -i "s/^queue [0-9]*/queue ${new_max}/" "$sub_file"
		if grep -q "^queue ${new_max}$" "$sub_file"; then
			log_status "  Updated $sub_file: queue $new_max — verified"
		else
			log_error "  FAILED to update $sub_file: 'queue $new_max' not found after sed"
			return 1
		fi
	fi

	(cd "$ANALYSIS_CODE_DIR" && \
	 git add "$(basename "$extras_file")" "$(basename "$sub_file")" 2>/dev/null && \
	 git commit -m "auto: update file_batch_max to $new_max for run_year $run_year_key (hadd chunking)" 2>/dev/null) && \
		log_status "  Git committed source changes" || \
		log_status "  Git commit skipped (no changes or not a repo)"
}

# ── Chunked hadd fallback ───────────────────────────────────────────────────
# When recursive_hadd fails, split input files into separate partN.root files.
# Updates record, source code, and .sub automatically.
# Sets FALLBACK_TOTAL_ENTRIES on success.
FALLBACK_TOTAL_ENTRIES=0
chunked_hadd_fallback() {
	local tid="$1" outds="$2" target_dir="$3" target_subdir="$4"
	local pre_count="$5" original_part="$6"
	shift 6
	local files=("$@")

	local n=${#files[@]}
	local num_chunks=$(( (n + HADD_CHUNK_MAX_FILES - 1) / HADD_CHUNK_MAX_FILES ))
	local file_prefix
	case "$target_subdir" in
		pp_2024)   file_prefix="data_pp24" ;;
		pbpb_2023) file_prefix="data_pbpb23" ;;
		pbpb_2024) file_prefix="data_pbpb24" ;;
		pbpb_2025) file_prefix="data_pbpb25" ;;
		*) log_error "task $tid: unknown target_subdir '$target_subdir'"; return 1 ;;
	esac

	local original_pnum
	original_pnum=$(echo "$original_part" | grep -oP '[0-9]+')

	log_status "task $tid: chunked fallback — splitting $n files into $num_chunks chunks of ≤$HADD_CHUNK_MAX_FILES"

	# Find current max part number on disk
	local max_existing=0
	for f in "$target_dir"/${file_prefix}_part*.root; do
		[[ -f "$f" ]] || continue
		local pnum
		pnum=$(basename "$f" | grep -oP 'part\K[0-9]+')
		if (( pnum > max_existing )); then max_existing=$pnum; fi
	done

	# Assign part numbers: first chunk keeps original_pnum, rest get max_existing+1, +2, ...
	local chunk_pnums=("$original_pnum")
	local next_new=$(( max_existing + 1 ))
	for (( i=1; i<num_chunks; i++ )); do
		if (( next_new == original_pnum )); then
			next_new=$(( next_new + 1 ))
		fi
		chunk_pnums+=("$next_new")
		next_new=$(( next_new + 1 ))
	done

	# hadd each chunk and validate
	local chunk_counts=()
	local total_entries=0
	for (( i=0; i<num_chunks; i++ )); do
		local start=$(( i * HADD_CHUNK_MAX_FILES ))
		local count=$HADD_CHUNK_MAX_FILES
		if (( start + count > n )); then count=$(( n - start )); fi
		local chunk=("${files[@]:$start:$count}")
		local pnum=${chunk_pnums[$i]}
		local chunk_out="${target_dir}/${file_prefix}_part${pnum}.root"

		log_status "task $tid: chunk $((i+1))/$num_chunks (${#chunk[@]} files) → ${file_prefix}_part${pnum}.root"

		rm -f "$chunk_out"
		if ! recursive_hadd "$chunk_out" "${chunk[@]}"; then
			log_error "task $tid: chunk $((i+1)) hadd failed"
			for (( j=0; j<=i; j++ )); do
				rm -f "${target_dir}/${file_prefix}_part${chunk_pnums[$j]}.root"
			done
			return 1
		fi

		local post_cnt pre_cnt
		post_cnt=$(count_entries_chain "$chunk_out")
		pre_cnt=$(count_entries_chain "${chunk[@]}")

		if [[ "$post_cnt" -ne "$pre_cnt" ]]; then
			log_error "task $tid: chunk $((i+1)) entry mismatch: pre=$pre_cnt post=$post_cnt"
			for (( j=0; j<=i; j++ )); do
				rm -f "${target_dir}/${file_prefix}_part${chunk_pnums[$j]}.root"
			done
			return 1
		fi

		chunk_counts+=("$post_cnt")
		total_entries=$(( total_entries + post_cnt ))
		log_status "task $tid: chunk $((i+1)) OK — $post_cnt entries"
	done

	# Validate total
	if [[ "$total_entries" -ne "$pre_count" ]]; then
		log_error "task $tid: TOTAL ENTRY MISMATCH: pre=$pre_count chunks_total=$total_entries"
		for (( j=0; j<num_chunks; j++ )); do
			rm -f "${target_dir}/${file_prefix}_part${chunk_pnums[$j]}.root"
		done
		return 1
	fi

	# Record all chunks
	for (( i=0; i<num_chunks; i++ )); do
		local pnum=${chunk_pnums[$i]}
		echo "${file_prefix}_part${pnum}.root | ${outds} | $(date +%F) | ${chunk_counts[$i]} entries (chunk $((i+1))/${num_chunks} of ${original_part})" >> "$RECORD_FILE"
	done
	log_status "task $tid: recorded $num_chunks chunks in $RECORD_FILE"

	# Update new max part number across all existing + newly created
	local new_max=$max_existing
	for pn in "${chunk_pnums[@]}"; do
		if (( pn > new_max )); then new_max=$pn; fi
	done

	# Update source code and .sub
	local code_info
	code_info=$(get_code_update_info "$target_subdir")
	if [[ -n "$code_info" ]]; then
		local extras_file sub_file run_year_key
		read -r extras_file sub_file run_year_key <<< "$code_info"
		update_source_for_new_max "$extras_file" "$sub_file" "$run_year_key" "$new_max"
	fi

	FALLBACK_TOTAL_ENTRIES=$total_entries
	log_status "task $tid: chunked fallback complete — $num_chunks chunks, max part=$new_max, total=$total_entries entries"
	return 0
}

# ── Download + hadd + validate for one task ──────────────────────────────────
# Split into two phases so retries only repeat hadd, never re-download.
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

	# ── Phase 1: Download (runs once; rucio skips already-downloaded files) ──
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

	# ── Phase 2: hadd (try single-file merge, fall back to separate chunks) ──
	log_status "task $tid: hadd (${#roots[@]} files) → $output_name.root"
	rm -f "$output_file"

	if recursive_hadd "$output_file" "${roots[@]}"; then
		# Single-file merge succeeded
		local post_count
		post_count=$(count_entries_chain "$output_file")
		log_status "task $tid: post-merge entries = $post_count"

		if [[ "$post_count" -ne "$pre_count" ]]; then
			log_error "task $tid: ENTRY MISMATCH for $output_name: pre=$pre_count post=$post_count"
			echo "${output_name}.root | ${outds} | $(date +%F) | ${post_count} entries !!! MISMATCH: source=${pre_count} merged=${post_count}" >> "$RECORD_FILE"
			return 1
		fi

		echo "${output_name}.root | ${outds} | $(date +%F) | ${post_count} entries" >> "$RECORD_FILE"
		log_status "task $tid: validated OK — ${post_count} entries."
	else
		# Single-file merge failed — chunked fallback
		rm -f "$output_file"
		log_status "task $tid: recursive hadd failed, using chunked fallback..."

		local part_str
		part_str=$(echo "$output_name" | grep -oP 'part[0-9]+')

		if ! chunked_hadd_fallback "$tid" "$outds" "$target_dir" "$target_subdir" "$pre_count" "$part_str" "${roots[@]}"; then
			return 1
		fi

		if [[ "$FALLBACK_TOTAL_ENTRIES" -ne "$pre_count" ]]; then
			log_error "task $tid: TOTAL ENTRY MISMATCH after chunked fallback: pre=$pre_count total=$FALLBACK_TOTAL_ENTRIES"
			return 1
		fi
		log_status "task $tid: chunked fallback OK — $FALLBACK_TOTAL_ENTRIES entries across multiple parts."
	fi

	log_status "task $tid: Cleaning up download dir."
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

# ── Parallel coordination ───────────────────────────────────────────────────

# Poll BigPanDA for pending tasks, mark ready ones.
refresh_task_states() {
	# Collect pending tasks from our list (brief lock)
	lock_state
	local pending_tids=()
	for tid in "${TASK_IDS[@]}"; do
		local state
		state=$(grep "^${tid} " "$STATE_FILE" | head -1 | awk '{print $2}')
		if [[ "$state" == "pending" ]]; then
			pending_tids+=("$tid")
		fi
	done
	unlock_state

	if [[ ${#pending_tids[@]} -gt 0 ]]; then
		log_status "Polling BigPanDA for ${#pending_tids[@]} pending task(s)..."

		for tid in "${pending_tids[@]}"; do
			local result
			result=$(query_task "$tid")
			IFS='|' read -r status pct taskname errmsg <<< "$result"

			case "$status" in
				done)
					log_status "  Task $tid: grid DONE (100%). Marking ready."
					set_task_state_if "$tid" "pending" "ready"
					;;
				finished)
					if (( $(echo "$pct >= $MIN_SUCCESS_PCT" | bc -l) )); then
						log_status "  Task $tid: grid FINISHED (${pct}%). Marking ready."
						set_task_state_if "$tid" "pending" "ready"
					else
						log_error "task $tid: FINISHED with low success (${pct}% < ${MIN_SUCCESS_PCT}%)"
						set_task_state_if "$tid" "pending" "failed"
					fi
					;;
				broken|aborted|failed|exhausted)
					log_error "task $tid: terminal grid status '$status'. $errmsg"
					set_task_state_if "$tid" "pending" "failed"
					;;
				api_error)
					log_status "  Task $tid: API error ($errmsg), will retry."
					;;
				*)
					log_status "  Task $tid: $status (${pct}%). Still running."
					;;
			esac
		done
	fi

	reclaim_stale_downloads
}

# Reset downloads that have been stuck longer than STALE_TIMEOUT
reclaim_stale_downloads() {
	lock_state
	local now=$(date +%s)
	local stale_tids=()
	local stale_hosts=()
	local stale_ages=()
	while IFS=' ' read -r tid state host epoch_ts; do
		if [[ "$state" == "downloading" && -n "$epoch_ts" ]]; then
			local age=$(( now - epoch_ts ))
			if (( age > STALE_TIMEOUT )); then
				stale_tids+=("$tid")
				stale_hosts+=("$host")
				stale_ages+=("$age")
			fi
		fi
	done < "$STATE_FILE"
	for i in "${!stale_tids[@]}"; do
		sed -i "s/^${stale_tids[$i]} downloading.*/${stale_tids[$i]} ready/" "$STATE_FILE"
	done
	unlock_state
	for i in "${!stale_tids[@]}"; do
		log_status "Reclaimed stale task ${stale_tids[$i]} (was ${stale_hosts[$i]}, ${stale_ages[$i]}s ago)"
	done
}

# Atomically claim one ready task from our task list. Echoes task ID or empty.
claim_one_ready_task() {
	lock_state
	local tid=""
	for t in "${TASK_IDS[@]}"; do
		local state
		state=$(grep "^${t} " "$STATE_FILE" | head -1 | awk '{print $2}')
		if [[ "$state" == "ready" ]]; then
			tid="$t"
			break
		fi
	done
	if [[ -n "$tid" ]]; then
		local now=$(date +%s)
		sed -i "s/^${tid} ready.*/${tid} downloading ${MY_HOSTNAME} ${now}/" "$STATE_FILE"
	fi
	unlock_state
	echo "$tid"
}

# Check if all tasks in our list are resolved (completed or failed)
all_tasks_resolved() {
	lock_state
	local unresolved=0
	for tid in "${TASK_IDS[@]}"; do
		local state
		state=$(grep "^${tid} " "$STATE_FILE" | head -1 | awk '{print $2}')
		case "$state" in
			pending|ready|downloading) unresolved=$((unresolved + 1)) ;;
		esac
	done
	unlock_state
	(( unresolved == 0 ))
}

# Print summary of task states
print_summary() {
	local completed=0 failed=0 pending=0 ready=0 downloading=0
	lock_state
	for tid in "${TASK_IDS[@]}"; do
		local state
		state=$(grep "^${tid} " "$STATE_FILE" | head -1 | awk '{print $2}')
		case "$state" in
			completed)   completed=$((completed + 1)) ;;
			failed)      failed=$((failed + 1)) ;;
			pending)     pending=$((pending + 1)) ;;
			ready)       ready=$((ready + 1)) ;;
			downloading) downloading=$((downloading + 1)) ;;
		esac
	done
	unlock_state
	log_status "Summary: ${completed} completed, ${failed} failed, ${downloading} downloading, ${ready} ready, ${pending} pending"
}

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

# Environment setup (ROOT for hadd, rucio for download)
log_status "Setting up environment..."
source ~/setup.sh || { echo "ERROR: setup.sh failed" >&2; exit 1; }
lsetup rucio      || { echo "ERROR: lsetup rucio failed" >&2; exit 1; }
log_status "Environment ready. root=$(command -v root) rucio=$(command -v rucio)"

init_state

log_status "════════════════════════════════════════════════════════════════"
log_status "Worker started: ${#TASK_IDS[@]} tasks, poll every ${POLL_INTERVAL_MIN}min"
log_status "  Status log: $STATUS_LOG"
log_status "  Error log:  $ERROR_LOG"
log_status "  State file: $STATE_FILE"
log_status "════════════════════════════════════════════════════════════════"

# Print initial state
lock_state
for tid in "${TASK_IDS[@]}"; do
	local_state=$(grep "^${tid} " "$STATE_FILE" | head -1 | awk '{print $2}')
	log_status "  $tid  ${OUTDS_MAP[$tid]:-<will resolve>}  [$local_state]"
done
unlock_state

# ── Worker loop ──────────────────────────────────────────────────────────────
while true; do
	# Step 1: Try to claim a ready task (fast, no API calls)
	claimed_tid=$(claim_one_ready_task)

	if [[ -z "$claimed_tid" ]]; then
		# Step 2: Nothing ready. Poll BigPanDA for newly completed grid tasks.
		refresh_task_states

		# Try again after refresh
		claimed_tid=$(claim_one_ready_task)
	fi

	if [[ -n "$claimed_tid" ]]; then
		# Check proxy before downloading
		if ! check_proxy; then
			log_error "Proxy expired. Releasing task $claimed_tid."
			set_task_state_if "$claimed_tid" "downloading" "ready"
			log_status "Sleeping 5min waiting for proxy renewal..."
			sleep 300
			continue
		fi

		# Resolve outDS
		outds="${OUTDS_MAP[$claimed_tid]:-}"
		if [[ -z "$outds" ]]; then
			outds=$(resolve_outds "$claimed_tid" "")
			if [[ -n "$outds" ]]; then
				OUTDS_MAP["$claimed_tid"]="$outds"
			fi
		fi

		if [[ -z "$outds" ]]; then
			log_error "task $claimed_tid: could not resolve outDS. Marking failed."
			set_task_state_if "$claimed_tid" "downloading" "failed"
			continue
		fi

		log_status "Claimed task $claimed_tid -> $outds"

		if process_task "$claimed_tid" "$outds"; then
			set_task_state_if "$claimed_tid" "downloading" "completed"
			log_status "Task $claimed_tid: COMPLETED successfully."
		else
			set_task_state_if "$claimed_tid" "downloading" "failed"
			log_error "task $claimed_tid: FAILED."
		fi

		print_summary
		continue  # Immediately check for more work
	fi

	# Step 3: Nothing ready even after polling. Check if all done.
	if all_tasks_resolved; then
		log_status ""
		log_status "════════════════════════════════════════════════════════════════"
		print_summary
		log_status "All tasks resolved. Worker exiting."
		log_status "════════════════════════════════════════════════════════════════"

		# Reorganize merging record via Claude Code (first worker to finish does this)
		if [[ -f "$RECORD_FILE" ]] && command -v claude &>/dev/null; then
			reorg_lock="${STATE_FILE}.reorg_lock"
			if (set -C; echo $$ > "$reorg_lock") 2>/dev/null; then
				log_status "Reorganizing $RECORD_FILE ..."
				RECORD_BACKUP="${RECORD_FILE}.bak"
				cp "$RECORD_FILE" "$RECORD_BACKUP"
				claude_err_file=$(mktemp)
				claude -p "$(cat << CLAUDEOF
You are reorganizing a bookkeeping log file for merged ROOT datasets.

File path: ${RECORD_FILE}

Each non-empty, non-comment line has the format:
  <output_file> | <source_subdir> | <date> | <N> entries[  !!! MISMATCH: source=<X> merged=<Y>]

The dataset key is derived from the output filename:
- For sequential-named files (contain _part): the prefix before _part (e.g. data_pbpb24, data_pp23).
- For custom-named files (no _part): the full basename without .root.

Task:
1. Read the file.
2. Group lines by dataset key. Discard existing comment lines (starting with #).
3. Within each group, sort by part number ascending (if present); otherwise by date ascending.
4. Between groups, sort alphabetically by dataset key.
5. Write the file back with: a blank line before each group (except the first),
   a comment header line  # <dataset_key>  before each group, then the sorted data lines.
6. Do not alter any data line's content, only reorder lines and add/replace headers.
CLAUDEOF
)" 2>"$claude_err_file"
				claude_exit=$?
				if [[ $claude_exit -eq 0 ]]; then
					log_status "  Record reorganized successfully"
					rm -f "$RECORD_BACKUP"
				else
					log_status "  Claude reorganization failed (exit $claude_exit) — record preserved as-is"
					cp "$RECORD_BACKUP" "$RECORD_FILE"
					rm -f "$RECORD_BACKUP"
				fi
				rm -f "$claude_err_file" "$reorg_lock"
			fi
		fi

		break
	fi

	# Step 4: Some tasks still pending or being downloaded by other workers.
	log_status "No tasks ready. Sleeping ${POLL_INTERVAL_MIN}min..."
	sleep "$POLL_INTERVAL"
done
