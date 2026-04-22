#!/usr/bin/env bash
# =============================================================================
# test_add_sample.sh
#
# Tests the incremental-sample use case with full-size STAR data:
#   Phase 1 — Full 4-step decoupled pipeline on 3 samples
#             (Div_100_S91, Div_101_S92, PARDOS_1_S1) to build a universe
#             and produce a 3-sample expression matrix.
#   Phase 2 — Add PARDOS_2_S2 WITHOUT re-running SCAN2 for the original 3.
#             Only SCAN1 and SCAN2 are run for the new sample; FINALIZE is
#             re-run with all 4 samples using the same universe, producing
#             a 4-sample expression matrix.
#
# Usage:
#   bash scripts/test_add_sample.sh [options]
#
#   --threads N      threads for each per-sample stage (default: 8)
#   --output-dir D   output directory (default: DATA_DIR/add_sample_test)
#   --intron         run with intron mode (-It 1)
#   --keep           keep output directory after the run
#
# Requirements:
#   - CIRI3_decoupled.jar at the repository root (run: bash scripts/build_jar.sh)
#   - java on PATH or in $CONDA_PREFIX/bin
#   - STAR output directories for all 4 samples under DATA_DIR
#   - Reference FASTA and GTF at the paths configured below
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
THREADS=8
KEEP=0
OUT_ROOT=""
INTRON_FLAG=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)    THREADS="$2";  shift 2 ;;
        --output-dir) OUT_ROOT="$2"; shift 2 ;;
        --intron)     INTRON_FLAG=(-It 1); shift ;;
        --keep)       KEEP=1;        shift   ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Paths - edit to match your environment
# ---------------------------------------------------------------------------
DATA_DIR=/pastel/tools/circRNA_tools/test_data
REF_FA=${DATA_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa
GTF_FILE=${DATA_DIR}/gencode.v32.primary_assembly.annotation.gtf

# First 3 samples form the initial universe; PARDOS_2_S2 is added in Phase 2.
INITIAL_SAMPLES=(
    "Div_100_S91"
    "Div_101_S92"
    "PARDOS_1_S1"
)
NEW_SAMPLE="PARDOS_2_S2"

# ---------------------------------------------------------------------------
# Derived paths
# ---------------------------------------------------------------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DECOUPLED_JAR="${REPO_ROOT}/CIRI3_decoupled.jar"

[[ -z "$OUT_ROOT" ]] && OUT_ROOT="${DATA_DIR}/add_sample_test"

SCAN1_DIR="${OUT_ROOT}/scan1"
UNIVERSE_DIR="${OUT_ROOT}/universe"
SCAN2_DIR="${OUT_ROOT}/scan2"
FINALIZE_3_DIR="${OUT_ROOT}/finalize_3samples"
FINALIZE_4_DIR="${OUT_ROOT}/finalize_4samples"
BENCH_DIR="${OUT_ROOT}/bench"

SCAN1_META_TSV="${UNIVERSE_DIR}/samples_scan1.tsv"
FINALIZE_3_TSV="${FINALIZE_3_DIR}/finalize_samples.tsv"
FINALIZE_4_TSV="${FINALIZE_4_DIR}/finalize_samples.tsv"

# shellcheck source=/dev/null
source "${REPO_ROOT}/scripts/_bench.sh"

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
CIRI3=("${JAVA_BIN}" -jar "${DECOUPLED_JAR}")

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
        rm -rf "$OUT_ROOT"
        info "Output directory removed."
    else
        info "Outputs kept at: $OUT_ROOT"
    fi
}
trap cleanup EXIT

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
info "=== Pre-flight checks ==="
[[ -f "$DECOUPLED_JAR" ]] || die "JAR not found: $DECOUPLED_JAR — run: bash scripts/build_jar.sh"
[[ -f "$REF_FA"        ]] || die "Reference FASTA not found: $REF_FA"
[[ -f "$GTF_FILE"      ]] || die "GTF not found: $GTF_FILE"
[[ -x "$JAVA_BIN"      ]] || die "Java not found (is the CIRI3 conda env active?)"

ALL_SAMPLES=("${INITIAL_SAMPLES[@]}" "$NEW_SAMPLE")
for S in "${ALL_SAMPLES[@]}"; do
    STAR_DIR="${DATA_DIR}/STAR_output_${S}"
    [[ -f "${STAR_DIR}/Chimeric.out.junction" ]] || die "Missing: ${STAR_DIR}/Chimeric.out.junction"
    [[ -f "${STAR_DIR}/Aligned.out.sam"       ]] || die "Missing: ${STAR_DIR}/Aligned.out.sam"
    [[ -f "${STAR_DIR}/bwa.sam"               ]] || die "Missing: ${STAR_DIR}/bwa.sam"
done
info "All input files found for ${#ALL_SAMPLES[@]} samples."
info "Output directory: ${OUT_ROOT}"
info "Decoupled JAR:    ${DECOUPLED_JAR}"

mkdir -p "$SCAN1_DIR" "$UNIVERSE_DIR" "$SCAN2_DIR" \
         "$FINALIZE_3_DIR" "$FINALIZE_4_DIR" "$BENCH_DIR"

> "$SCAN1_META_TSV"
> "$FINALIZE_3_TSV"

# ============================================================================
# PHASE 1: Full 4-step pipeline with the initial 3 samples
# ============================================================================
echo ""
echo "############################################################"
echo "#  PHASE 1: Build universe from 3 samples (full pipeline)  #"
echo "############################################################"

# --- Stage 1: SCAN1 (Div_100_S91, Div_101_S92, PARDOS_1_S1) ---
info "--- Phase 1 / Stage 1: SCAN1 ---"
SCAN1_IDX=0
for SAMPLE_ID in "${INITIAL_SAMPLES[@]}"; do
    SCAN1_IDX=$((SCAN1_IDX+1))
    STAR_DIR="${DATA_DIR}/STAR_output_${SAMPLE_ID}"
    TRIPLE="${STAR_DIR}/Chimeric.out.junction,${STAR_DIR}/Aligned.out.sam,${STAR_DIR}/bwa.sam"
    BWA_SAM="${STAR_DIR}/bwa.sam"
    mkdir -p "${SCAN1_DIR}/${SAMPLE_ID}"
    OUT_PREFIX="${SCAN1_DIR}/${SAMPLE_ID}/${SAMPLE_ID}"
    META="${OUT_PREFIX}.scan1_meta"

    if [[ -s "$META" ]]; then
        info "  [SKIP] SCAN1 already done for ${SAMPLE_ID}"
    else
        rm -f "${OUT_PREFIX}BSJ"* "${BWA_SAM}BSJ"*
        info "  SCAN1: ${SAMPLE_ID}"
        scan1_stdout="${BENCH_DIR}/scan1_${SAMPLE_ID}.stdout"
        bench_run "$(printf '10_scan1_%02d_%s' "${SCAN1_IDX}" "${SAMPLE_ID}")" \
            "${CIRI3[@]}" SCAN1 \
                -I "${TRIPLE}" \
                -O "${OUT_PREFIX}" \
                -F "${REF_FA}" \
                -A "${GTF_FILE}" \
                -T "${THREADS}" -Ma 1 -S 2 "${INTRON_FLAG[@]}" \
            2>&1 | tee "${scan1_stdout}" \
            | grep -iE "SCAN1|scan.*completed|meta|time|Mapped|Exception|Error|BrokenBarrier" || true
        if [[ ! -s "$META" ]]; then
            echo "[DEBUG] Last 30 lines of SCAN1 stdout (${scan1_stdout}):"
            tail -30 "${scan1_stdout}" 2>/dev/null || echo "  (no stdout captured)"
            java_log="${OUT_PREFIX}.log"
            [[ -s "$java_log" ]] && { echo "[DEBUG] Java log:"; tail -30 "${java_log}"; }
        fi
    fi

    check_exists "SCAN1 meta (${SAMPLE_ID})" "${META}"

    if [[ ! -s "$META" ]]; then
        SPLIT_NUM=0
    else
        SPLIT_NUM=$(grep "^fileSplitNum=" "${META}" | cut -d= -f2)
        local_fail=0
        for bsj_i in $(seq 1 "${SPLIT_NUM}"); do
            [[ -f "${OUT_PREFIX}BSJ${bsj_i}" ]] || { fail "Missing BSJ file: ${OUT_PREFIX}BSJ${bsj_i}"; local_fail=1; }
        done
        [[ $local_fail -eq 0 ]] && ok "SCAN1 BSJ files present for ${SAMPLE_ID} (${SPLIT_NUM} splits)"
    fi

    echo -e "${BWA_SAM}\t${META}" >> "$SCAN1_META_TSV"
done

# --- Stage 2: BUILD_UNIVERSE (3 samples) ---
info "--- Phase 1 / Stage 2: BUILD_UNIVERSE ---"
UNIVERSE_FILE="${UNIVERSE_DIR}/cohort.universe"
if [[ -s "$UNIVERSE_FILE" ]]; then
    info "  [SKIP] Universe already exists"
else
    bench_run "20_build_universe_3samples" \
        "${CIRI3[@]}" BUILD_UNIVERSE \
            -I "${SCAN1_META_TSV}" \
            -F "${REF_FA}" \
            -O "${UNIVERSE_DIR}/cohort" \
        2>&1 | grep -E "Universe|circRNA|time" || true
fi
check_exists "Universe file (3 samples)" "${UNIVERSE_FILE}"
UNIVERSE_CIRCS=$(grep -c "^chr" "${UNIVERSE_FILE}" || true)
info "Universe (3 samples): ${UNIVERSE_CIRCS} circRNA candidates."

# --- Stage 3: SCAN2 (Div_100_S91, Div_101_S92, PARDOS_1_S1) ---
info "--- Phase 1 / Stage 3: SCAN2 (initial 3 samples) ---"
SCAN2_IDX=0
for SAMPLE_ID in "${INITIAL_SAMPLES[@]}"; do
    SCAN2_IDX=$((SCAN2_IDX+1))
    STAR_DIR="${DATA_DIR}/STAR_output_${SAMPLE_ID}"
    TRIPLE="${STAR_DIR}/Chimeric.out.junction,${STAR_DIR}/Aligned.out.sam,${STAR_DIR}/bwa.sam"
    BWA_SAM="${STAR_DIR}/bwa.sam"
    SCAN1_PREFIX="${SCAN1_DIR}/${SAMPLE_ID}/${SAMPLE_ID}"
    META="${SCAN1_PREFIX}.scan1_meta"
    SPLIT_NUM=$(grep "^fileSplitNum=" "${META}" | cut -d= -f2)
    mkdir -p "${SCAN2_DIR}/${SAMPLE_ID}"
    OUT_PREFIX="${SCAN2_DIR}/${SAMPLE_ID}/${SAMPLE_ID}"
    FSJ_COUNTS="${OUT_PREFIX}.fsj_counts"

    if [[ -s "$FSJ_COUNTS" ]]; then
        info "  [SKIP] SCAN2 already done for ${SAMPLE_ID}"
    else
        info "  SCAN2: ${SAMPLE_ID}"
        bench_run "$(printf '30_scan2_%02d_%s' "${SCAN2_IDX}" "${SAMPLE_ID}")" \
            "${CIRI3[@]}" SCAN2 \
                -I "${TRIPLE}" \
                -CU "${UNIVERSE_FILE}" \
                -SM "${META}" \
                -O "${OUT_PREFIX}" \
                -F "${REF_FA}" \
                -T "${THREADS}" -Ma 1 "${INTRON_FLAG[@]}" \
            2>&1 | grep -E "scan|FSJ|BSJ|time" || true
    fi
    check_exists "SCAN2 FSJ counts (${SAMPLE_ID})" "${FSJ_COUNTS}"
    echo -e "${BWA_SAM}\t${FSJ_COUNTS}\t${SPLIT_NUM}\t${SAMPLE_ID}\t${SCAN1_PREFIX}" >> "$FINALIZE_3_TSV"
done

# --- Stage 4: FINALIZE (3 samples) ---
info "--- Phase 1 / Stage 4: FINALIZE (3 samples) ---"
FINAL_3_BSJ="${FINALIZE_3_DIR}/result_3samples.BSJ_Matrix"
if [[ -s "$FINAL_3_BSJ" ]]; then
    info "  [SKIP] 3-sample FINALIZE already done"
else
    bench_run "40_finalize_3samples" \
        "${CIRI3[@]}" FINALIZE \
            -I "${FINALIZE_3_TSV}" \
            -CU "${UNIVERSE_FILE}" \
            -F "${REF_FA}" \
            -O "${FINALIZE_3_DIR}/result_3samples" \
            -A "${GTF_FILE}" \
            -S 2 "${INTRON_FLAG[@]}" \
        2>&1 | grep -E "FINALIZE|Summary|Matrix|circRNA|time" || true
fi
check_exists "3-sample BSJ_Matrix" "${FINAL_3_BSJ}"
check_exists "3-sample FSJ_Matrix" "${FINALIZE_3_DIR}/result_3samples.FSJ_Matrix"
CIRCS_3=$(tail -n +2 "${FINAL_3_BSJ}" | wc -l)
info "3-sample result: $(matrix_dims "${FINAL_3_BSJ}")"

# ============================================================================
# PHASE 2: Add PARDOS_2_S2 — SCAN1 + SCAN2 only, reuse existing universe
# ============================================================================
echo ""
echo "###########################################################"
echo "#  PHASE 2: Add PARDOS_2_S2 — skip SCAN2 for original 3  #"
echo "###########################################################"

info "New sample: ${NEW_SAMPLE}"
info "Universe reused from Phase 1 — SCAN2 is NOT re-run for the original 3 samples."

NEW_STAR_DIR="${DATA_DIR}/STAR_output_${NEW_SAMPLE}"
NEW_TRIPLE="${NEW_STAR_DIR}/Chimeric.out.junction,${NEW_STAR_DIR}/Aligned.out.sam,${NEW_STAR_DIR}/bwa.sam"
NEW_BWA_SAM="${NEW_STAR_DIR}/bwa.sam"
NEW_SCAN1_PREFIX="${SCAN1_DIR}/${NEW_SAMPLE}/${NEW_SAMPLE}"
NEW_META="${NEW_SCAN1_PREFIX}.scan1_meta"

# --- Phase 2 / Stage 1: SCAN1 for PARDOS_2_S2 ---
info "--- Phase 2 / Stage 1: SCAN1 (${NEW_SAMPLE}) ---"
mkdir -p "${SCAN1_DIR}/${NEW_SAMPLE}"
if [[ -s "$NEW_META" ]]; then
    info "  [SKIP] SCAN1 already done for ${NEW_SAMPLE}"
else
    rm -f "${NEW_SCAN1_PREFIX}BSJ"* "${NEW_BWA_SAM}BSJ"*
    info "  SCAN1: ${NEW_SAMPLE}"
    scan1_stdout="${BENCH_DIR}/scan1_${NEW_SAMPLE}.stdout"
    bench_run "50_scan1_${NEW_SAMPLE}" \
        "${CIRI3[@]}" SCAN1 \
            -I "${NEW_TRIPLE}" \
            -O "${NEW_SCAN1_PREFIX}" \
            -F "${REF_FA}" \
            -A "${GTF_FILE}" \
            -T "${THREADS}" -Ma 1 -S 2 "${INTRON_FLAG[@]}" \
        2>&1 | tee "${scan1_stdout}" \
        | grep -iE "SCAN1|scan.*completed|meta|time|Mapped|Exception|Error|BrokenBarrier" || true
    if [[ ! -s "$NEW_META" ]]; then
        echo "[DEBUG] Last 30 lines of SCAN1 stdout (${scan1_stdout}):"
        tail -30 "${scan1_stdout}" 2>/dev/null || echo "  (no stdout captured)"
        java_log="${NEW_SCAN1_PREFIX}.log"
        [[ -s "$java_log" ]] && { echo "[DEBUG] Java log:"; tail -30 "${java_log}"; }
    fi
fi

check_exists "SCAN1 meta (${NEW_SAMPLE})" "${NEW_META}"

if [[ ! -s "$NEW_META" ]]; then
    NEW_SPLIT_NUM=0
else
    NEW_SPLIT_NUM=$(grep "^fileSplitNum=" "${NEW_META}" | cut -d= -f2)
    new_local_fail=0
    for bsj_i in $(seq 1 "${NEW_SPLIT_NUM}"); do
        [[ -f "${NEW_SCAN1_PREFIX}BSJ${bsj_i}" ]] || { fail "Missing BSJ file: ${NEW_SCAN1_PREFIX}BSJ${bsj_i}"; new_local_fail=1; }
    done
    [[ $new_local_fail -eq 0 ]] && ok "SCAN1 BSJ files present for ${NEW_SAMPLE} (${NEW_SPLIT_NUM} splits)"
fi

# --- Phase 2 / Stage 2: SCAN2 for PARDOS_2_S2 against the existing universe ---
info "--- Phase 2 / Stage 2: SCAN2 (${NEW_SAMPLE}) against existing 3-sample universe ---"
mkdir -p "${SCAN2_DIR}/${NEW_SAMPLE}"
NEW_SCAN2_PREFIX="${SCAN2_DIR}/${NEW_SAMPLE}/${NEW_SAMPLE}"
NEW_FSJ_COUNTS="${NEW_SCAN2_PREFIX}.fsj_counts"

if [[ -s "$NEW_FSJ_COUNTS" ]]; then
    info "  [SKIP] SCAN2 already done for ${NEW_SAMPLE}"
else
    bench_run "60_scan2_${NEW_SAMPLE}" \
        "${CIRI3[@]}" SCAN2 \
            -I "${NEW_TRIPLE}" \
            -CU "${UNIVERSE_FILE}" \
            -SM "${NEW_META}" \
            -O "${NEW_SCAN2_PREFIX}" \
            -F "${REF_FA}" \
            -T "${THREADS}" -Ma 1 "${INTRON_FLAG[@]}" \
        2>&1 | grep -E "scan|FSJ|BSJ|time" || true
fi
check_exists "SCAN2 FSJ counts (${NEW_SAMPLE})" "${NEW_FSJ_COUNTS}"

# --- Phase 2 / Stage 3: FINALIZE with all 4 samples ---
info "--- Phase 2 / Stage 3: FINALIZE (all 4 samples) ---"
FINAL_4_BSJ="${FINALIZE_4_DIR}/result_4samples.BSJ_Matrix"

# Build 4-sample FINALIZE TSV: original 3 entries + PARDOS_2_S2
cp "$FINALIZE_3_TSV" "$FINALIZE_4_TSV"
echo -e "${NEW_BWA_SAM}\t${NEW_FSJ_COUNTS}\t${NEW_SPLIT_NUM}\t${NEW_SAMPLE}\t${NEW_SCAN1_PREFIX}" \
    >> "$FINALIZE_4_TSV"

if [[ -s "$FINAL_4_BSJ" ]]; then
    info "  [SKIP] 4-sample FINALIZE already done"
else
    bench_run "70_finalize_4samples" \
        "${CIRI3[@]}" FINALIZE \
            -I "${FINALIZE_4_TSV}" \
            -CU "${UNIVERSE_FILE}" \
            -FM "${FINAL_3_BSJ}" \
            -F "${REF_FA}" \
            -O "${FINALIZE_4_DIR}/result_4samples" \
            -A "${GTF_FILE}" \
            -S 2 "${INTRON_FLAG[@]}" \
        2>&1 | grep -E "FINALIZE|Summary|Frozen|Matrix|circRNA|time" || true
fi
check_exists "4-sample BSJ_Matrix" "${FINAL_4_BSJ}"
check_exists "4-sample FSJ_Matrix" "${FINALIZE_4_DIR}/result_4samples.FSJ_Matrix"
CIRCS_4=$(tail -n +2 "${FINAL_4_BSJ}" | wc -l)
info "4-sample result: $(matrix_dims "${FINAL_4_BSJ}")"

# ============================================================================
# Verification
# ============================================================================
echo ""
info "=== Verification ==="

# --- Column structure ---
COL_3=$(head -1 "${FINAL_3_BSJ}" | awk '{print NF-1}')
[[ "$COL_3" -eq 3 ]] && ok "3-sample BSJ_Matrix has exactly 3 sample columns" \
    || fail "3-sample BSJ_Matrix: expected 3 columns, got $COL_3"

COL_4=$(head -1 "${FINAL_4_BSJ}" | awk '{print NF-1}')
[[ "$COL_4" -eq 4 ]] && ok "4-sample BSJ_Matrix has exactly 4 sample columns" \
    || fail "4-sample BSJ_Matrix: expected 4 columns, got $COL_4"

head -1 "${FINAL_4_BSJ}" | grep -q "${NEW_SAMPLE}" \
    && ok "${NEW_SAMPLE} column present in 4-sample BSJ_Matrix" \
    || fail "${NEW_SAMPLE} column MISSING from 4-sample BSJ_Matrix"

for S in "${INITIAL_SAMPLES[@]}"; do
    head -1 "${FINAL_4_BSJ}" | grep -q "$S" \
        && ok "Column '$S' retained in 4-sample matrix" \
        || fail "Column '$S' missing from 4-sample matrix"
done

# --- circRNA row identity ---
# The 4-sample matrix must have the same rows as the 3-sample matrix.
# Adding a new sample must not change which circRNAs are reported.
IDS_3=$(tail -n +2 "${FINAL_3_BSJ}" | awk '{print $1}' | sort)
IDS_4=$(tail -n +2 "${FINAL_4_BSJ}" | awk '{print $1}' | sort)

if [[ "$CIRCS_3" -eq "$CIRCS_4" ]]; then
    ok "circRNA row count identical: ${CIRCS_3} in both 3-sample and 4-sample matrices"
else
    fail "circRNA row count differs: 3-sample=${CIRCS_3}, 4-sample=${CIRCS_4} (expected equal)"
fi

if [[ "$IDS_3" == "$IDS_4" ]]; then
    ok "circRNA IDs identical between 3-sample and 4-sample matrices"
else
    ONLY_3=$(comm -23 <(echo "$IDS_3") <(echo "$IDS_4") | wc -l)
    ONLY_4=$(comm -13 <(echo "$IDS_3") <(echo "$IDS_4") | wc -l)
    fail "circRNA IDs differ: ${ONLY_3} only in 3-sample, ${ONLY_4} only in 4-sample"
    info "  (first 5 only in 3-sample):"; comm -23 <(echo "$IDS_3") <(echo "$IDS_4") | head -5 | sed 's/^/    /' || true
    info "  (first 5 only in 4-sample):"; comm -13 <(echo "$IDS_3") <(echo "$IDS_4") | head -5 | sed 's/^/    /' || true
fi

# --- Count identity for original 3 samples ---
# BSJ and FSJ counts for the original samples must be byte-for-byte identical.
# Only the PARDOS_2_S2 column should differ.
for MATRIX_SUFFIX in BSJ_Matrix FSJ_Matrix; do
    F3="${FINALIZE_3_DIR}/result_3samples.${MATRIX_SUFFIX}"
    F4="${FINALIZE_4_DIR}/result_4samples.${MATRIX_SUFFIX}"
    [[ -s "$F3" && -s "$F4" ]] || continue

    # Build sorted (circRNA_ID, col_value) pairs for each original sample
    HDR_3=$(head -1 "${F3}")
    HDR_4=$(head -1 "${F4}")
    mismatch=0
    for S in "${INITIAL_SAMPLES[@]}"; do
        # Column index in 3-sample matrix (1-based field; col 1 = circRNA_ID)
        C3=$(echo "$HDR_3" | tr '\t' '\n' | grep -n "^${S}$" | cut -d: -f1)
        C4=$(echo "$HDR_4" | tr '\t' '\n' | grep -n "^${S}$" | cut -d: -f1)
        if [[ -z "$C3" || -z "$C4" ]]; then
            fail "${MATRIX_SUFFIX}: column '${S}' not found in one of the matrices"
            mismatch=1; continue
        fi
        DATA_3=$(tail -n +2 "${F3}" | awk -v c="$C3" 'BEGIN{OFS="\t"}{print $1,$c}' | sort)
        DATA_4=$(tail -n +2 "${F4}" | awk -v c="$C4" 'BEGIN{OFS="\t"}{print $1,$c}' | sort)
        if [[ "$DATA_3" == "$DATA_4" ]]; then
            ok "${MATRIX_SUFFIX}: counts for '${S}' identical in 3-sample and 4-sample matrices"
        else
            NDIFF=$(diff <(echo "$DATA_3") <(echo "$DATA_4") | grep -c "^[<>]" || true)
            fail "${MATRIX_SUFFIX}: counts for '${S}' differ (${NDIFF} differing lines)"
            diff <(echo "$DATA_3") <(echo "$DATA_4") | head -10 || true
            mismatch=1
        fi
    done
done

# ============================================================================
# Benchmark report + summary
# ============================================================================
bench_report

echo ""
echo "========================================================"
echo "  3-sample matrix : ${FINAL_3_BSJ}"
echo "  4-sample matrix : ${FINAL_4_BSJ}"
echo "========================================================"
echo "  TEST SUMMARY"
echo "========================================================"
printf "  3-sample pipeline : %d circRNAs\n" "${CIRCS_3}"
printf "  4-sample pipeline : %d circRNAs\n" "${CIRCS_4}"
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
