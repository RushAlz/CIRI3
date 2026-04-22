#!/usr/bin/env bash
# =============================================================================
# test_add_sample.sh
#
# Tests the incremental-sample use case:
#   Phase 1 — Full 4-step pipeline with 3 samples (sample1, sample2, sample3)
#             to create a universe and produce a 3-sample expression matrix.
#   Phase 2 — Add PARDOS_2_S2 (sample4.sam) WITHOUT re-running SCAN2 for the
#             original 3 samples.  Only SCAN1 and SCAN2 are run for the new
#             sample; FINALIZE is then re-run with all 4 samples against the
#             same universe, producing a 4-sample expression matrix.
#
# Usage:
#   bash scripts/test_add_sample.sh [--threads N] [--keep]
#
#   --threads N   number of threads for per-sample stages (default: 1)
#   --keep        do not delete the working directory after the test
#
# Requirements:
#   - CIRI3_decoupled.jar at the repository root (build with scripts/build_jar.sh)
#   - java on PATH or in $CONDA_PREFIX/bin
#   - sample{1,2,3,4}.sam in data/circRNA/Mutiple/
#   - Reference FASTA at data/circRNA/ref.fa
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
THREADS=1
KEEP=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --keep)    KEEP=1; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/_bench.sh"

JAR="$REPO_ROOT/CIRI3_decoupled.jar"
DATA_DIR="$REPO_ROOT/data/circRNA/Mutiple"
REF_FA="$REPO_ROOT/data/circRNA/ref.fa"

WORK_DIR="$REPO_ROOT/test_add_sample_$$"
SCAN1_DIR="$WORK_DIR/scan1"
UNIVERSE_DIR="$WORK_DIR/universe"
SCAN2_DIR="$WORK_DIR/scan2"
FINALIZE_3_DIR="$WORK_DIR/finalize_3samples"
FINALIZE_4_DIR="$WORK_DIR/finalize_4samples"
BENCH_DIR="$WORK_DIR/bench"

PASS=0
FAIL=0

# ---------------------------------------------------------------------------
# Java binary
# ---------------------------------------------------------------------------
if [[ -n "${CONDA_PREFIX:-}" ]] && [[ -x "${CONDA_PREFIX}/bin/java" ]]; then
    JAVA_BIN="${CONDA_PREFIX}/bin/java"
else
    JAVA_BIN="$(command -v java)"
fi

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
info()  { echo "[INFO]  $*"; }
ok()    { echo "[PASS]  $*"; PASS=$((PASS+1)); }
fail()  { echo "[FAIL]  $*"; FAIL=$((FAIL+1)); }
die()   { echo "[ERROR] $*" >&2; exit 1; }

check_exists() {
    local label="$1" f="$2"
    if [[ -s "$f" ]]; then ok "$label exists and is non-empty"
    else fail "$label missing or empty: $f"; fi
}

matrix_dims() {
    local f="$1"
    local cols rows
    cols=$(head -1 "$f" | awk '{print NF-1}')
    rows=$(tail -n +2 "$f" | wc -l)
    echo "${rows} circRNAs x ${cols} samples"
}

cleanup() {
    if [[ $KEEP -eq 0 ]]; then
        rm -rf "$WORK_DIR"
        info "Working directory removed."
    else
        info "Working directory kept at: $WORK_DIR"
    fi
}
trap cleanup EXIT

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
info "=== Pre-flight checks ==="
[[ -f "$JAR"     ]] || die "JAR not found: $JAR — run: bash scripts/build_jar.sh"
[[ -d "$DATA_DIR" ]] || die "Test data not found at $DATA_DIR"
[[ -f "$REF_FA"  ]] || die "Reference FASTA not found at $REF_FA"
[[ -x "$JAVA_BIN"  ]] || die "Java not found (is the CIRI3 conda env active?)"

for s in 1 2 3 4; do
    [[ -f "$DATA_DIR/sample${s}.sam" ]] || die "Missing test SAM: $DATA_DIR/sample${s}.sam"
done
info "JAR: $JAR"
info "Found all 4 sample SAM files."

mkdir -p "$SCAN1_DIR" "$UNIVERSE_DIR" "$SCAN2_DIR" \
         "$FINALIZE_3_DIR" "$FINALIZE_4_DIR" "$BENCH_DIR"

# TSV files built incrementally
SCAN1_META_TSV="$UNIVERSE_DIR/samples_scan1.tsv"
FINALIZE_3_TSV="$FINALIZE_3_DIR/finalize_samples.tsv"
FINALIZE_4_TSV="$FINALIZE_4_DIR/finalize_samples.tsv"
> "$SCAN1_META_TSV"
> "$FINALIZE_3_TSV"

# Shorthand for running the JAR
CIRI3="${JAVA_BIN} -jar ${JAR}"

# ============================================================================
# PHASE 1: Full 4-step pipeline with the initial 3 samples
# ============================================================================
echo ""
echo "############################################################"
echo "#  PHASE 1: Build universe from 3 samples (full pipeline)  #"
echo "############################################################"

INITIAL_SAMPLES=("sample1" "sample2" "sample3")

# --- Stage 1: SCAN1 (sample1, sample2, sample3) ---
info "--- Phase 1 / Stage 1: SCAN1 ---"
for SAMPLE_NAME in "${INITIAL_SAMPLES[@]}"; do
    SAM="$DATA_DIR/${SAMPLE_NAME}.sam"
    mkdir -p "$SCAN1_DIR/${SAMPLE_NAME}"
    OUT_PREFIX="$SCAN1_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}"
    info "  SCAN1: $SAMPLE_NAME"
    bench_run "10_scan1_${SAMPLE_NAME}" \
        ${CIRI3} SCAN1 \
            -I "$SAM" \
            -O "$OUT_PREFIX" \
            -F "$REF_FA" \
            -T "$THREADS" \
            -S 2 \
        2>&1 | grep -E "scan|completed|meta|time|ERROR|Exception" || true
    check_exists "SCAN1 meta ($SAMPLE_NAME)" "${OUT_PREFIX}.scan1_meta"
    echo -e "$SAM\t${OUT_PREFIX}.scan1_meta" >> "$SCAN1_META_TSV"
done

# --- Stage 2: BUILD_UNIVERSE (3 samples) ---
info "--- Phase 1 / Stage 2: BUILD_UNIVERSE ---"
bench_run "20_build_universe_3samples" \
    ${CIRI3} BUILD_UNIVERSE \
        -I "$SCAN1_META_TSV" \
        -F "$REF_FA" \
        -O "$UNIVERSE_DIR/cohort" \
    2>&1 | grep -E "Universe|circRNA|time|ERROR|Exception" || true
check_exists "Universe file (3 samples)" "$UNIVERSE_DIR/cohort.universe"
UNIVERSE_LINES=$(grep -c "^chr" "$UNIVERSE_DIR/cohort.universe" || true)
info "Universe (3 samples) contains $UNIVERSE_LINES circRNA candidates."

# --- Stage 3: SCAN2 (sample1, sample2, sample3) ---
info "--- Phase 1 / Stage 3: SCAN2 (initial 3 samples) ---"
for SAMPLE_NAME in "${INITIAL_SAMPLES[@]}"; do
    SAM="$DATA_DIR/${SAMPLE_NAME}.sam"
    SCAN1_META="$SCAN1_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.scan1_meta"
    SPLIT_NUM=$(grep "^fileSplitNum=" "$SCAN1_META" | cut -d= -f2)
    SCAN1_PREFIX="$SCAN1_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}"
    mkdir -p "$SCAN2_DIR/${SAMPLE_NAME}"
    OUT_PREFIX="$SCAN2_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}"
    info "  SCAN2: $SAMPLE_NAME"
    bench_run "30_scan2_${SAMPLE_NAME}" \
        ${CIRI3} SCAN2 \
            -I "$SAM" \
            -CU "$UNIVERSE_DIR/cohort.universe" \
            -SM "$SCAN1_META" \
            -O "$OUT_PREFIX" \
            -F "$REF_FA" \
            -T "$THREADS" \
        2>&1 | grep -E "scan|completed|FSJ|time|ERROR|Exception" || true
    check_exists "SCAN2 FSJ counts ($SAMPLE_NAME)" "${OUT_PREFIX}.fsj_counts"
    echo -e "$SAM\t${OUT_PREFIX}.fsj_counts\t${SPLIT_NUM}\t${SAMPLE_NAME}\t${SCAN1_PREFIX}" \
        >> "$FINALIZE_3_TSV"
done

# --- Stage 4: FINALIZE (3 samples) ---
info "--- Phase 1 / Stage 4: FINALIZE (3 samples) ---"
bench_run "40_finalize_3samples" \
    ${CIRI3} FINALIZE \
        -I "$FINALIZE_3_TSV" \
        -CU "$UNIVERSE_DIR/cohort.universe" \
        -F "$REF_FA" \
        -O "$FINALIZE_3_DIR/result_3samples" \
        -S 2 \
    2>&1 | grep -E "FINALIZE|Summary|Matrix|circRNA|time|ERROR|Exception" || true
check_exists "3-sample BSJ_Matrix" "$FINALIZE_3_DIR/result_3samples.BSJ_Matrix"
check_exists "3-sample FSJ_Matrix" "$FINALIZE_3_DIR/result_3samples.FSJ_Matrix"
info "3-sample result: $(matrix_dims "$FINALIZE_3_DIR/result_3samples.BSJ_Matrix")"

# ============================================================================
# PHASE 2: Add PARDOS_2_S2 (sample4) — SCAN1 + SCAN2 only for the new sample
# ============================================================================
echo ""
echo "###########################################################"
echo "#  PHASE 2: Add PARDOS_2_S2 — skip SCAN2 for original 3  #"
echo "###########################################################"

NEW_SAMPLE_NAME="PARDOS_2_S2"
NEW_SAM="$DATA_DIR/sample4.sam"

info "New sample: $NEW_SAMPLE_NAME (data: sample4.sam)"
info "Universe reused from Phase 1 — no SCAN2 re-run for sample1/sample2/sample3."

# --- Stage 1: SCAN1 for PARDOS_2_S2 ---
info "--- Phase 2 / Stage 1: SCAN1 (PARDOS_2_S2) ---"
mkdir -p "$SCAN1_DIR/${NEW_SAMPLE_NAME}"
NEW_SCAN1_PREFIX="$SCAN1_DIR/${NEW_SAMPLE_NAME}/${NEW_SAMPLE_NAME}"
bench_run "50_scan1_${NEW_SAMPLE_NAME}" \
    ${CIRI3} SCAN1 \
        -I "$NEW_SAM" \
        -O "$NEW_SCAN1_PREFIX" \
        -F "$REF_FA" \
        -T "$THREADS" \
        -S 2 \
    2>&1 | grep -E "scan|completed|meta|time|ERROR|Exception" || true
check_exists "SCAN1 meta ($NEW_SAMPLE_NAME)" "${NEW_SCAN1_PREFIX}.scan1_meta"

# --- Stage 2: SCAN2 for PARDOS_2_S2 against the existing 3-sample universe ---
info "--- Phase 2 / Stage 2: SCAN2 (PARDOS_2_S2) against existing universe ---"
NEW_SCAN1_META="${NEW_SCAN1_PREFIX}.scan1_meta"
NEW_SPLIT_NUM=$(grep "^fileSplitNum=" "$NEW_SCAN1_META" | cut -d= -f2)
mkdir -p "$SCAN2_DIR/${NEW_SAMPLE_NAME}"
NEW_SCAN2_PREFIX="$SCAN2_DIR/${NEW_SAMPLE_NAME}/${NEW_SAMPLE_NAME}"
bench_run "60_scan2_${NEW_SAMPLE_NAME}" \
    ${CIRI3} SCAN2 \
        -I "$NEW_SAM" \
        -CU "$UNIVERSE_DIR/cohort.universe" \
        -SM "$NEW_SCAN1_META" \
        -O "$NEW_SCAN2_PREFIX" \
        -F "$REF_FA" \
        -T "$THREADS" \
    2>&1 | grep -E "scan|completed|FSJ|time|ERROR|Exception" || true
check_exists "SCAN2 FSJ counts ($NEW_SAMPLE_NAME)" "${NEW_SCAN2_PREFIX}.fsj_counts"

# --- Stage 3: FINALIZE with all 4 samples ---
info "--- Phase 2 / Stage 3: FINALIZE (all 4 samples) ---"

# Build the 4-sample FINALIZE TSV: reuse original 3 entries + add PARDOS_2_S2
cp "$FINALIZE_3_TSV" "$FINALIZE_4_TSV"
echo -e "$NEW_SAM\t${NEW_SCAN2_PREFIX}.fsj_counts\t${NEW_SPLIT_NUM}\t${NEW_SAMPLE_NAME}\t${NEW_SCAN1_PREFIX}" \
    >> "$FINALIZE_4_TSV"

bench_run "70_finalize_4samples" \
    ${CIRI3} FINALIZE \
        -I "$FINALIZE_4_TSV" \
        -CU "$UNIVERSE_DIR/cohort.universe" \
        -F "$REF_FA" \
        -O "$FINALIZE_4_DIR/result_4samples" \
        -S 2 \
    2>&1 | grep -E "FINALIZE|Summary|Matrix|circRNA|time|ERROR|Exception" || true
check_exists "4-sample BSJ_Matrix" "$FINALIZE_4_DIR/result_4samples.BSJ_Matrix"
check_exists "4-sample FSJ_Matrix" "$FINALIZE_4_DIR/result_4samples.FSJ_Matrix"
info "4-sample result: $(matrix_dims "$FINALIZE_4_DIR/result_4samples.BSJ_Matrix")"

# ============================================================================
# Verification
# ============================================================================
echo ""
info "=== Verification ==="

# 3-sample matrix must have exactly 3 sample columns
COL_3=$(head -1 "$FINALIZE_3_DIR/result_3samples.BSJ_Matrix" | awk '{print NF-1}')
if [[ "$COL_3" -eq 3 ]]; then
    ok "3-sample BSJ_Matrix has exactly 3 sample columns"
else
    fail "3-sample BSJ_Matrix: expected 3 sample columns, got $COL_3"
fi

# 4-sample matrix must have exactly 4 sample columns
COL_4=$(head -1 "$FINALIZE_4_DIR/result_4samples.BSJ_Matrix" | awk '{print NF-1}')
if [[ "$COL_4" -eq 4 ]]; then
    ok "4-sample BSJ_Matrix has exactly 4 sample columns"
else
    fail "4-sample BSJ_Matrix: expected 4 sample columns, got $COL_4"
fi

# PARDOS_2_S2 column must appear in the 4-sample matrix header
if head -1 "$FINALIZE_4_DIR/result_4samples.BSJ_Matrix" | grep -q "PARDOS_2_S2"; then
    ok "PARDOS_2_S2 column is present in the 4-sample BSJ_Matrix"
else
    fail "PARDOS_2_S2 column is MISSING from the 4-sample BSJ_Matrix"
fi

# 4-sample matrix must have at least as many circRNAs as the 3-sample matrix
ROWS_3=$(tail -n +2 "$FINALIZE_3_DIR/result_3samples.BSJ_Matrix" | wc -l)
ROWS_4=$(tail -n +2 "$FINALIZE_4_DIR/result_4samples.BSJ_Matrix" | wc -l)
if [[ "$ROWS_4" -ge "$ROWS_3" ]]; then
    ok "4-sample matrix circRNA count ($ROWS_4) >= 3-sample count ($ROWS_3)"
else
    fail "4-sample matrix circRNA count ($ROWS_4) < 3-sample count ($ROWS_3)"
fi

info "Universe file used: $UNIVERSE_DIR/cohort.universe ($UNIVERSE_LINES candidates)"

# Original 3 sample columns must still be present in the 4-sample matrix
for S in sample1 sample2 sample3; do
    if head -1 "$FINALIZE_4_DIR/result_4samples.BSJ_Matrix" | grep -q "$S"; then
        ok "Column '$S' retained in 4-sample matrix"
    else
        fail "Column '$S' missing from 4-sample matrix"
    fi
done

# ============================================================================
# Benchmark report + summary
# ============================================================================
bench_report

echo ""
echo "========================================================"
echo "  3-sample matrix : $FINALIZE_3_DIR/result_3samples.BSJ_Matrix"
echo "  4-sample matrix : $FINALIZE_4_DIR/result_4samples.BSJ_Matrix"
echo "========================================================"
echo "  TEST SUMMARY"
echo "========================================================"
echo "  PASSED : $PASS"
echo "  FAILED : $FAIL"
echo "========================================================"

if [[ $FAIL -eq 0 ]]; then
    echo "  ALL TESTS PASSED"
    exit 0
else
    echo "  SOME TESTS FAILED — see output above"
    exit 1
fi
