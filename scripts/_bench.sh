# shellcheck shell=bash
# =============================================================================
# scripts/_bench.sh
#
# Lightweight benchmarking helpers for the CIRI3 test scripts.
#
# Provides two functions:
#   bench_run  LABEL CMD [ARGS...]   - run a command; record wall time,
#                                       CPU time, %CPU, and peak RSS.
#   bench_report                     - print a table of all recorded stages
#                                       with a total row.
#
# Requires `BENCH_DIR` to be set (a writable directory) BEFORE sourcing or
# before calling bench_run. Uses GNU `time` (/usr/bin/time) when available
# for accurate peak-RSS measurement; falls back to bash `SECONDS` for
# wall-clock only when GNU time is missing.
#
# On error exit codes, bench_run propagates the command's exit status via
# PIPESTATUS so callers that pipe through `grep` / `tee` can still get it.
# =============================================================================

# --- Locate GNU time --------------------------------------------------------
if [[ -x /usr/bin/time ]]; then
    BENCH_GTIME=/usr/bin/time
elif command -v gtime >/dev/null 2>&1; then
    BENCH_GTIME="$(command -v gtime)"
else
    BENCH_GTIME=""
fi

# --- Format used when GNU time is available ---------------------------------
#   wall_sec / user_sec / sys_sec : seconds, floating point
#   max_rss_kb                    : KiB
#   cpu_pct                       : percentage string (e.g. "647%")
#   exit                          : exit code
BENCH_FMT='wall_sec=%e user_sec=%U sys_sec=%S max_rss_kb=%M cpu_pct=%P exit=%x'

bench_run() {
    local label="$1"; shift
    if [[ -z "${BENCH_DIR:-}" ]]; then
        echo "[bench] BENCH_DIR not set; running without measurement" >&2
        "$@"
        return $?
    fi
    mkdir -p "${BENCH_DIR}"
    local bench_file="${BENCH_DIR}/${label}.bench"

    if [[ -n "${BENCH_GTIME}" ]]; then
        "${BENCH_GTIME}" -f "${BENCH_FMT}" -o "${bench_file}" "$@"
        return $?
    fi

    # Fallback: wall-clock only via bash SECONDS.
    local _bench_start=${SECONDS}
    local _bench_status=0
    "$@" || _bench_status=$?
    local _bench_elapsed=$(( SECONDS - _bench_start ))
    printf 'wall_sec=%d user_sec=- sys_sec=- max_rss_kb=- cpu_pct=- exit=%d\n' \
        "${_bench_elapsed}" "${_bench_status}" > "${bench_file}"
    return ${_bench_status}
}

# Format a number of seconds as Hh Mm Ss when >= 60s, else "Ns".
_bench_fmt_time() {
    awk -v s="$1" 'BEGIN {
        if (s == "-") { print "-"; exit }
        s = s + 0
        if (s < 60) { printf "%.1fs", s; exit }
        h = int(s / 3600); m = int((s - h*3600) / 60); r = s - h*3600 - m*60
        if (h > 0) printf "%dh%02dm%02ds", h, m, r
        else       printf "%dm%02ds",       m,    r
    }'
}

_bench_fmt_gb() {
    awk -v k="$1" 'BEGIN {
        if (k == "-") { print "-"; exit }
        printf "%.2f", k / 1024 / 1024
    }'
}

bench_report() {
    echo ""
    echo "================================================================"
    echo "  BENCHMARK: wall time & peak RAM per stage"
    echo "================================================================"
    if [[ -z "${BENCH_DIR:-}" ]] || ! compgen -G "${BENCH_DIR}/*.bench" >/dev/null 2>&1; then
        echo "  (no benchmark data)"
        echo "================================================================"
        return
    fi
    if [[ -z "${BENCH_GTIME}" ]]; then
        echo "  NOTE: GNU time not found; only wall-clock captured."
        echo "        Install 'time' (apt install time / conda install time)"
        echo "        for peak-RSS measurement."
    fi
    printf "  %-32s %12s %12s %12s %10s\n" "Stage" "Wall" "CPU" "CPU%" "MaxRSS(GB)"
    printf "  %-32s %12s %12s %12s %10s\n" "--------------------------------" "------------" "------------" "------------" "----------"

    local total_wall=0
    local peak_rss=0
    # Stable order: sort benchmark files by name (prefixed with zero-padded
    # sequence numbers by the caller).
    local bench_file
    while IFS= read -r bench_file; do
        local label wall_sec user_sec sys_sec max_rss_kb cpu_pct exit
        label=$(basename "$bench_file" .bench)
        # shellcheck disable=SC1090
        wall_sec="-" user_sec="-" sys_sec="-" max_rss_kb="-" cpu_pct="-" exit="-"
        eval "$(cat "$bench_file")"
        local cpu_total
        cpu_total=$(awk -v u="$user_sec" -v s="$sys_sec" 'BEGIN {
            if (u == "-" || s == "-") { print "-"; exit }
            printf "%.1f", u + s
        }')
        printf "  %-32s %12s %12s %12s %10s\n" \
            "$label" \
            "$(_bench_fmt_time "$wall_sec")" \
            "$(_bench_fmt_time "$cpu_total")" \
            "${cpu_pct}" \
            "$(_bench_fmt_gb "$max_rss_kb")"
        if [[ "$wall_sec" != "-" ]]; then
            total_wall=$(awk -v a="$total_wall" -v b="$wall_sec" 'BEGIN { printf "%.1f", a + b }')
        fi
        if [[ "$max_rss_kb" != "-" ]] && (( max_rss_kb > peak_rss )); then
            peak_rss=$max_rss_kb
        fi
    done < <(ls -1 "${BENCH_DIR}"/*.bench | sort)

    printf "  %-32s %12s %12s %12s %10s\n" "--------------------------------" "------------" "------------" "------------" "----------"
    local rss_total
    if (( peak_rss > 0 )); then
        rss_total=$(_bench_fmt_gb "$peak_rss")
    else
        rss_total="-"
    fi
    printf "  %-32s %12s %12s %12s %10s\n" \
        "TOTAL (sum wall / peak RSS)" \
        "$(_bench_fmt_time "$total_wall")" \
        "-" "-" \
        "$rss_total"
    echo "================================================================"
}
