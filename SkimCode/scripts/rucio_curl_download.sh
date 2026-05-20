#!/bin/bash
# rucio_curl_download.sh — Download rucio datasets using REST API + curl.
# Workaround for broken CVMFS Python (BNL Stratum 1 returning 503 for
# sft.cern.ch Python 3.11 stdlib). Uses zero CVMFS Python — only curl,
# system python3 for JSON parsing, and grid CA certs from /etc.
#
# Usage:
#   ./rucio_curl_download.sh may2026_skim.txt
#   ./rucio_curl_download.sh user.yuhang.TrigRates.dimuon.pp2024data.May2026.v1.part3._EXT0

set -o pipefail

DATA_BASE="/usatlas/u/yuhanguo/usatlasdata/dimuon_data"
PROXY_CERT="${X509_USER_PROXY:-/tmp/x509up_u$(id -u)}"
CA_OPT="--capath /etc/grid-security/certificates/"
AUTH_HOST="https://voatlasrucio-auth-prod.cern.ch"
RUCIO_HOST="https://voatlasrucio-server-prod.cern.ch"
RUCIO_ACCOUNT="yuhang"
MAX_RETRIES=3
PARALLEL_DL=4

ts() { date "+%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] $*"; }
err() { echo "[$(ts)] ERROR: $*" >&2; }

# ── Parse input ──────────────────────────────────────────────────────────────
DATASETS=()
while [[ $# -gt 0 ]]; do
    if [[ -f "$1" ]]; then
        while IFS= read -r line; do
            [[ "$line" =~ ^[[:space:]]*# ]] && continue
            [[ -z "$line" ]] && continue
            outds=$(echo "$line" | grep -oE 'user\.yuhang\.[A-Za-z0-9._]+_EXT0' | head -1)
            [[ -n "$outds" ]] && DATASETS+=("$outds")
        done < "$1"
    else
        DATASETS+=("$1")
    fi
    shift
done

if [[ ${#DATASETS[@]} -eq 0 ]]; then
    echo "Usage: $0 FILE_OR_DATASETS..." >&2; exit 1
fi

# ── Auth ─────────────────────────────────────────────────────────────────────
get_token() {
    curl -s -i $CA_OPT --connect-timeout 10 --max-time 20 \
        --cert "$PROXY_CERT" --key "$PROXY_CERT" \
        -H "X-Rucio-Account: $RUCIO_ACCOUNT" \
        "${AUTH_HOST}/auth/x509_proxy" 2>/dev/null \
    | grep "X-Rucio-Auth-Token:" | sed 's/.*: //' | tr -d '\r'
}

log "Authenticating..."
TOKEN=$(get_token)
if [[ -z "$TOKEN" ]]; then
    err "Failed to get rucio auth token. Check proxy: voms-proxy-info -all"
    exit 1
fi
log "Auth OK (account: $RUCIO_ACCOUNT)"

# ── outDS → target dir/name mapping ─────────────────────────────────────────
map_outds() {
    local outds="$1"
    local part dir prefix
    part=$(echo "$outds" | grep -oP 'part\d+')
    if   [[ "$outds" == *PbPb2023* ]]; then dir="pbpb_2023"; prefix="data_pbpb23"
    elif [[ "$outds" == *PbPb2024* ]]; then dir="pbpb_2024"; prefix="data_pbpb24"
    elif [[ "$outds" == *PbPb2025* ]]; then dir="pbpb_2025"; prefix="data_pbpb25"
    elif [[ "$outds" == *pp2024*   ]]; then dir="pp_2024";   prefix="data_pp24"
    else echo "UNKNOWN UNKNOWN"; return 1; fi
    echo "${dir} ${prefix}_${part}"
}

# ── Get replicas for a dataset ───────────────────────────────────────────────
get_replicas() {
    local outds="$1"
    curl -s $CA_OPT --connect-timeout 10 --max-time 60 \
        --cert "$PROXY_CERT" --key "$PROXY_CERT" \
        -H "X-Rucio-Auth-Token: $TOKEN" \
        -H "Content-Type: application/json" \
        -X POST \
        -d "{\"dids\": [{\"scope\": \"user.yuhang\", \"name\": \"$outds\"}]}" \
        "${RUCIO_HOST}/replicas/list" 2>/dev/null
}

# ── Pick ranked download URLs for each file ──────────────────────────────────
# Output: name\tbytes\tadler32\turl1|url2|url3 (pipe-separated, ranked)
pick_urls() {
    /usr/bin/python3 -c "
import json, sys

def to_https(url):
    if url.startswith('davs://'):
        return url.replace('davs://', 'https://')
    if url.startswith('root://'):
        import re
        return re.sub(r'^root://([^:]+):\d+/', r'https://\1/', url)
    return url

for line in sys.stdin:
    line = line.strip()
    if not line: continue
    r = json.loads(line)
    name = r['name']
    pfns = list(r.get('pfns', {}).keys())
    ranked = []
    # Tier 1: CERN EOS
    for u in pfns:
        if 'eosatlas.cern.ch' in u:
            ranked.append(to_https(u))
    # Tier 2: other davs/https sites (most reliable for curl)
    for u in pfns:
        if u.startswith('davs://') and 'eosatlas.cern.ch' not in u:
            ranked.append(to_https(u))
    # Tier 3: root:// converted to https
    for u in pfns:
        if u.startswith('root://') and 'eosatlas.cern.ch' not in u:
            ranked.append(to_https(u))
    # Deduplicate preserving order
    seen = set(); uniq = []
    for u in ranked:
        if u not in seen: seen.add(u); uniq.append(u)
    print(f'{name}\t{r[\"bytes\"]}\t{r.get(\"adler32\",\"\")}\t{\"|\".join(uniq)}')
"
}

# ── Download one file, trying multiple URLs ──────────────────────────────────
download_file() {
    local outpath="$1" expected_size="$2"
    shift 2
    local urls=("$@")
    for url in "${urls[@]}"; do
        local host
        host=$(echo "$url" | sed 's|https://||;s|/.*||')
        curl -s -L $CA_OPT --connect-timeout 15 --max-time 600 \
            --cert "$PROXY_CERT" --key "$PROXY_CERT" \
            -o "$outpath" "$url" 2>/dev/null
        local actual_size
        actual_size=$(stat -c%s "$outpath" 2>/dev/null || echo 0)
        if [[ "$actual_size" -eq "$expected_size" ]]; then
            return 0
        fi
        log "    $host failed (got $actual_size bytes), trying next replica..."
        rm -f "$outpath"
    done
    return 1
}

# ── Process one dataset ──────────────────────────────────────────────────────
process_dataset() {
    local outds="$1"
    local mapping target_subdir output_name
    mapping=$(map_outds "$outds")
    if [[ "$mapping" == "UNKNOWN UNKNOWN" ]]; then
        err "cannot map $outds"; return 1
    fi
    target_subdir=$(echo "$mapping" | cut -d' ' -f1)
    output_name=$(echo "$mapping" | cut -d' ' -f2)
    local dl_dir="${DATA_BASE}/${target_subdir}/${outds}"
    mkdir -p "$dl_dir"

    log "[$outds] Fetching replica list..."
    local replica_info
    replica_info=$(get_replicas "$outds" | pick_urls)
    if [[ -z "$replica_info" ]]; then
        TOKEN=$(get_token)
        replica_info=$(get_replicas "$outds" | pick_urls)
    fi
    if [[ -z "$replica_info" ]]; then
        err "[$outds] No replicas found"; return 1
    fi

    local nfiles
    nfiles=$(echo "$replica_info" | wc -l)
    log "[$outds] $nfiles files to download → $target_subdir/"

    local failed=0 downloaded=0
    while IFS=$'\t' read -r fname fsize fadler furls_joined; do
        local fpath="${dl_dir}/${fname}"
        if [[ -f "$fpath" ]] && [[ "$(stat -c%s "$fpath" 2>/dev/null)" == "$fsize" ]]; then
            log "  $fname — already downloaded, skipping"
            ((downloaded++))
            continue
        fi
        # Split pipe-separated URLs into array
        IFS='|' read -ra url_array <<< "$furls_joined"
        local primary_host
        primary_host=$(echo "${url_array[0]}" | sed 's|https://||;s|/.*||')
        log "  $fname ($(numfmt --to=iec "$fsize")) ← $primary_host (${#url_array[@]} replicas)"
        if download_file "$fpath" "$fsize" "${url_array[@]}"; then
            ((downloaded++))
        else
            err "  $fname — download failed from all ${#url_array[@]} replicas"
            ((failed++))
        fi
    done <<< "$replica_info"

    if [[ $failed -gt 0 ]]; then
        err "[$outds] $failed/$nfiles files failed"
        return 1
    fi
    log "[$outds] All $downloaded files downloaded to $dl_dir/"
    return 0
}

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
log "═══════════════════════════════════════════════════════════"
log "rucio_curl_download: ${#DATASETS[@]} datasets"
log "═══════════════════════════════════════════════════════════"

completed=0 failed=0
for ds in "${DATASETS[@]}"; do
    if process_dataset "$ds"; then
        ((completed++))
    else
        ((failed++))
    fi
done

log "═══════════════════════════════════════════════════════════"
log "Done: $completed completed, $failed failed out of ${#DATASETS[@]}"
log "═══════════════════════════════════════════════════════════"

if [[ $completed -gt 0 ]]; then
    log ""
    log "Files are in per-dataset subdirectories. To merge (once ROOT is available):"
    log "  source ~/setup.sh && hadd -f <output>.root <input_dir>/*.root"
fi

exit $((failed > 0 ? 1 : 0))
