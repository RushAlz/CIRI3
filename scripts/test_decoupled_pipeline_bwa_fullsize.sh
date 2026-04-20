#!/usr/bin/env bash
# =============================================================================
# test_decoupled_pipeline_bwa_fullsize.sh
#
# Runs the original joint (-W 1) BWA-only CIRI3 pipeline AND the four-stage
# decoupled pipeline (SCAN1 -> BUILD_UNIVERSE -> SCAN2 -> FINALIZE) on
# full-size BWA data, then verifies that both pipelines produce identical
# BSJ and FSJ matrices.
#
# The original (joint) pipeline runs the stock CIRI3_Java_1.8.0.jar so the
# comparison is against the published ground truth. The decoupled pipeline
# runs CIRI3_decoupled.jar, the jar built from this repo's src/ tree.
#
# Each stage is skipped if its outputs already exist, so the script is safe
# to re-run after a partial failure.
#
# Usage:
#   bash scripts/test_decoupled_pipeline_bwa_fullsize.sh [options]
#
#   --threads N          threads for each per-sample stage (default: 8)
#   --output-dir D       output directory (default: DATA_DIR/decoupled_bwa_comparison)
#   --intron             run with intron mode (-It 1) in both pipelines
#   --use-current-joint  run the joint pipeline from CIRI3_decoupled.jar (this
#                        repo's current source) instead of CIRI3_Java_1.8.0.jar.
#                        Use this to test the decoupled decomposition against
#                        the exact same BSJ-detection code, isolating any
#                        version-drift differences in the published jar.
#   --keep               keep output directory after the run
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
THREADS=8
KEEP=0
OUT_ROOT=""
INTRON_FLAG=()
USE_CURRENT_JOINT=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)            THREADS="$2";  shift 2 ;;
        --output-dir)         OUT_ROOT="$2"; shift 2 ;;
        --intron)             INTRON_FLAG=(-It 1); shift ;;
        --use-current-joint)  USE_CURRENT_JOINT=1; shift ;;
        --keep)               KEEP=1;        shift   ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Paths - edit to match your environment
# ---------------------------------------------------------------------------
DATA_DIR=/pastel/tools/circRNA_tools/test_data
REF_FA=${DATA_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa
GTF_FILE=${DATA_DIR}/gencode.v32.primary_assembly.annotation.gtf

SAMPLES=(
    "Div_100_S91"
    "Div_101_S92"
    "PARDOS_1_S1"
    "PARDOS_2_S2"
)

# ---------------------------------------------------------------------------
# Derived paths
# ---------------------------------------------------------------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

ORIGINAL_JAR="${REPO_ROOT}/CIRI3_Java_1.8.0.jar"       # published ground truth
DECOUPLED_JAR="${REPO_ROOT}/CIRI3_decoupled.jar"       # built from this repo

[[ -z "$OUT_ROOT" ]] && OUT_ROOT="${DATA_DIR}/decoupled_bwa_comparison"

ORIG_DIR="${OUT_ROOT}/original"
DECOUPLED_DIR="${OUT_ROOT}/decoupled"
SCAN1_DIR="${DECOUPLED_DIR}/scan1"
UNIVERSE_DIR="${DECOUPLED_DIR}/universe"
SCAN2_DIR="${DECOUPLED_DIR}/scan2"
FINALIZE_DIR="${DECOUPLED_DIR}/finalize"

SCAN1_META_TSV="${UNIVERSE_DIR}/samples_scan1.tsv"
FINALIZE_TSV="${FINALIZE_DIR}/finalize_samples.tsv"
BENCH_DIR="${OUT_ROOT}/bench"

# shellcheck source=/dev/null
source "${REPO_ROOT}/scripts/_bench.sh"

PASS=0
FAIL=0

# ---------------------------------------------------------------------------
# Java invocations - prefer the conda-managed JVM when available
# ---------------------------------------------------------------------------
if [[ -n "${CONDA_PREFIX:-}" ]] && [[ -x "${CONDA_PREFIX}/bin/java" ]]; then
    JAVA_BIN="${CONDA_PREFIX}/bin/java"
else
    JAVA_BIN="$(command -v java)"
fi

JAVA_NEW=(${JAVA_BIN} -jar "${DECOUPLED_JAR}")
if [[ $USE_CURRENT_JOINT -eq 1 ]]; then
    JAVA_ORIG=(${JAVA_BIN} -jar "${DECOUPLED_JAR}")
    JOINT_LABEL="CIRI3_decoupled.jar (current source, -W 1)"
else
    JAVA_ORIG=(${JAVA_BIN} -jar "${ORIGINAL_JAR}")
    JOINT_LABEL="CIRI3_Java_1.8.0.jar (published, -W 1)"
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
# 0. Pre-flight checks
# ---------------------------------------------------------------------------
info "=== Pre-flight checks ==="
[[ -f "$REF_FA"        ]] || die "Reference FASTA not found: $REF_FA"
[[ -f "$GTF_FILE"      ]] || die "GTF not found: $GTF_FILE"
[[ -x "$JAVA_BIN"      ]] || die "Java not found at: $JAVA_BIN (is the CIRI3 conda env active?)"
[[ -f "$ORIGINAL_JAR"  ]] || die "Original jar not found: $ORIGINAL_JAR"
[[ -f "$DECOUPLED_JAR" ]] || die "Decoupled jar not found: $DECOUPLED_JAR (run: scripts/build_jar.sh or see README)"

for S in "${SAMPLES[@]}"; do
    [[ -f "${DATA_DIR}/${S}.sam" ]] || die "Missing: ${DATA_DIR}/${S}.sam"
done
info "All input files found for ${#SAMPLES[@]} samples."
info "Output directory: ${OUT_ROOT}"
info "Joint pipeline: ${JOINT_LABEL}"
info "Decoupled jar:  ${DECOUPLED_JAR}"

mkdir -p "$ORIG_DIR" "$SCAN1_DIR" "$UNIVERSE_DIR" "$SCAN2_DIR" "$FINALIZE_DIR" "$BENCH_DIR"

# ---------------------------------------------------------------------------
# 1. ORIGINAL pipeline (-W 1, BWA-only) from CIRI3_Java_1.8.0.jar
# ---------------------------------------------------------------------------
info "=== Stage 0: ORIGINAL pipeline (CIRI3_Java_1.8.0.jar, -W 1, BWA) ==="

ORIG_BSJ="${ORIG_DIR}/result.BSJ_Matrix"
if [[ -s "$ORIG_BSJ" ]]; then
    info "  [SKIP] Original results exist at ${ORIG_DIR}/result"
else
    ORIG_TSV="${ORIG_DIR}/samples.tsv"
    > "$ORIG_TSV"
    for S in "${SAMPLES[@]}"; do
        echo "${DATA_DIR}/${S}.sam" >> "$ORIG_TSV"
    done

    bench_run "00_original_joint" \
        "${JAVA_ORIG[@]}" \
            -I "${ORIG_TSV}" \
            -O "${ORIG_DIR}/result" \
            -F "${REF_FA}" \
            -A "${GTF_FILE}" \
            -W 1 -T "${THREADS}" -S 0 "${INTRON_FLAG[@]}" \
        2>&1 | tee "${ORIG_DIR}/run.log" \
        | grep -E "CIRI3|scan|completed|circRNA|Mapped|time|Exception|Error|^\t?at |DIAG" || true
fi

if [[ ! -s "${ORIG_DIR}/result.BSJ_Matrix" ]]; then
    echo "[ERROR] Original pipeline produced no BSJ_Matrix — aborting." >&2
    echo "        Inspect: ${ORIG_DIR}/run.log" >&2
    tail -40 "${ORIG_DIR}/run.log" >&2 || true
    exit 1
fi

check_exists "Original BSJ_Matrix" "${ORIG_DIR}/result.BSJ_Matrix"
check_exists "Original FSJ_Matrix" "${ORIG_DIR}/result.FSJ_Matrix"
ORIG_CIRCS=$(tail -n +2 "${ORIG_DIR}/result.BSJ_Matrix" | wc -l)
info "Original pipeline: ${ORIG_CIRCS} circRNAs."

# ---------------------------------------------------------------------------
# 2. DECOUPLED pipeline (CIRI3_decoupled.jar, BWA-only)
# ---------------------------------------------------------------------------
info "=== Decoupled pipeline (CIRI3_decoupled.jar, BWA) ==="

> "$SCAN1_META_TSV"
> "$FINALIZE_TSV"

# --- Stage 1: SCAN1 (per sample) ---
info "--- Stage 1: SCAN1 ---"
SCAN1_IDX=0
for S in "${SAMPLES[@]}"; do
    SCAN1_IDX=$((SCAN1_IDX+1))
    BWA_SAM="${DATA_DIR}/${S}.sam"
    OUT_PREFIX="${SCAN1_DIR}/${S}"
    META="${OUT_PREFIX}.scan1_meta"

    if [[ -s "$META" ]]; then
        info "  [SKIP] SCAN1 already done for $S"
    else
        # Wipe any stale BSJ files from a previous run (different -T) so that
        # SCAN2's on-disk BSJ enumeration only sees files from this run.
        rm -f "${BWA_SAM}BSJ"*
        info "  SCAN1: $S"
        bench_run "$(printf '10_scan1_%02d_%s' "${SCAN1_IDX}" "${S}")" \
            "${JAVA_NEW[@]}" SCAN1 \
                -I "${BWA_SAM}" \
                -O "${OUT_PREFIX}" \
                -F "${REF_FA}" \
                -A "${GTF_FILE}" \
                -T "${THREADS}" -S 0 "${INTRON_FLAG[@]}" \
            2>&1 | grep -E "scan|meta|time|Mapped" || true
    fi

    check_exists "SCAN1 meta ($S)" "${META}"

    SPLIT_NUM=$(grep "^fileSplitNum=" "${META}" | cut -d= -f2)
    local_fail=0
    for i in $(seq 1 "${SPLIT_NUM}"); do
        [[ -f "${BWA_SAM}BSJ${i}" ]] || { fail "Missing BSJ file: ${BWA_SAM}BSJ${i}"; local_fail=1; }
    done
    [[ $local_fail -eq 0 ]] && ok "SCAN1 BSJ files present for $S (${SPLIT_NUM} splits)"

    echo -e "${BWA_SAM}\t${META}" >> "$SCAN1_META_TSV"
done

# --- Stage 2: BUILD_UNIVERSE ---
info "--- Stage 2: BUILD_UNIVERSE ---"
UNIVERSE_FILE="${UNIVERSE_DIR}/cohort.universe"
if [[ -s "$UNIVERSE_FILE" ]]; then
    info "  [SKIP] Universe already exists"
else
    bench_run "20_build_universe" \
        "${JAVA_NEW[@]}" BUILD_UNIVERSE \
            -I "${SCAN1_META_TSV}" \
            -F "${REF_FA}" \
            -O "${UNIVERSE_DIR}/cohort" \
        2>&1 | grep -E "Universe|circRNA|time" || true
fi
check_exists "Universe file" "${UNIVERSE_FILE}"
UNIVERSE_CIRCS=$(grep -c "^chr" "${UNIVERSE_FILE}" || true)
info "Universe: ${UNIVERSE_CIRCS} circRNA candidates."

# --- Stage 3: SCAN2 (per sample) ---
info "--- Stage 3: SCAN2 ---"
SCAN2_IDX=0
for S in "${SAMPLES[@]}"; do
    SCAN2_IDX=$((SCAN2_IDX+1))
    BWA_SAM="${DATA_DIR}/${S}.sam"
    META="${SCAN1_DIR}/${S}.scan1_meta"
    SPLIT_NUM=$(grep "^fileSplitNum=" "${META}" | cut -d= -f2)
    OUT_PREFIX="${SCAN2_DIR}/${S}"
    FSJ_COUNTS="${OUT_PREFIX}.fsj_counts"

    if [[ -s "$FSJ_COUNTS" ]]; then
        info "  [SKIP] SCAN2 already done for $S"
    else
        info "  SCAN2: $S"
        bench_run "$(printf '30_scan2_%02d_%s' "${SCAN2_IDX}" "${S}")" \
            "${JAVA_NEW[@]}" SCAN2 \
                -I "${BWA_SAM}" \
                -CU "${UNIVERSE_FILE}" \
                -O "${OUT_PREFIX}" \
                -F "${REF_FA}" \
                -T "${THREADS}" "${INTRON_FLAG[@]}" \
            2>&1 | grep -E "scan|FSJ|BSJ|time" || true
    fi
    check_exists "SCAN2 FSJ counts ($S)" "${FSJ_COUNTS}"

    echo -e "${BWA_SAM}\t${FSJ_COUNTS}\t${SPLIT_NUM}\t${S}" >> "$FINALIZE_TSV"
done

# --- Stage 4: FINALIZE ---
info "--- Stage 4: FINALIZE ---"
FINAL_BSJ="${FINALIZE_DIR}/result.BSJ_Matrix"
if [[ -s "$FINAL_BSJ" ]]; then
    info "  [SKIP] FINALIZE already done"
else
    bench_run "40_finalize" \
        "${JAVA_NEW[@]}" FINALIZE \
            -I "${FINALIZE_TSV}" \
            -CU "${UNIVERSE_FILE}" \
            -F "${REF_FA}" \
            -O "${FINALIZE_DIR}/result" \
            -A "${GTF_FILE}" \
            -S 0 "${INTRON_FLAG[@]}" \
        2>&1 | grep -E "FINALIZE|Summary|Matrix|circRNA|time" || true
fi
check_exists "Decoupled BSJ_Matrix" "${FINAL_BSJ}"
check_exists "Decoupled FSJ_Matrix" "${FINALIZE_DIR}/result.FSJ_Matrix"
DECOUPLED_CIRCS=$(tail -n +2 "${FINAL_BSJ}" | wc -l)
info "Decoupled pipeline: ${DECOUPLED_CIRCS} circRNAs."

# ---------------------------------------------------------------------------
# 3. Compare outputs
# ---------------------------------------------------------------------------
info "=== Comparing outputs ==="

normalise_matrix() {
    local f="$1"
    awk 'NR==1{
        printf "circRNA_ID"
        for(i=2;i<=NF;i++) printf "\ts%d", i-1
        printf "\n"; next
    }
    { print | "sort" }' "$f"
}

A_BSJ=$(normalise_matrix "${ORIG_DIR}/result.BSJ_Matrix")
B_BSJ=$(normalise_matrix "${FINAL_BSJ}")
if [[ "$A_BSJ" == "$B_BSJ" ]]; then
    ok "BSJ_Matrix: original and decoupled are IDENTICAL (${ORIG_CIRCS} circRNAs)"
else
    fail "BSJ_Matrix: original and decoupled DIFFER"
    { diff <(echo "$A_BSJ") <(echo "$B_BSJ") || true; } | head -40 || true
fi

A_FSJ=$(normalise_matrix "${ORIG_DIR}/result.FSJ_Matrix")
B_FSJ=$(normalise_matrix "${FINALIZE_DIR}/result.FSJ_Matrix")
if [[ "$A_FSJ" == "$B_FSJ" ]]; then
    ok "FSJ_Matrix: original and decoupled are IDENTICAL (${ORIG_CIRCS} circRNAs)"
else
    fail "FSJ_Matrix: original and decoupled DIFFER"
    { diff <(echo "$A_FSJ") <(echo "$B_FSJ") || true; } | head -40 || true
fi

info "Checking universe coverage..."
ORIG_IDS=$(tail -n +2 "${ORIG_DIR}/result.BSJ_Matrix" | awk '{print $1}' | sort)
UNIV_IDS=$(grep "^chr" "${UNIVERSE_FILE}" | awk '{printf "%s:%s|%s\n", $1, $2, $3}' | sort)
MISSING=$(comm -23 <(echo "$ORIG_IDS") <(echo "$UNIV_IDS") | wc -l)
if [[ "$MISSING" -eq 0 ]]; then
    ok "Universe coverage: all original circRNAs present in universe"
else
    fail "Universe coverage: $MISSING original circRNAs missing from universe"
    { comm -23 <(echo "$ORIG_IDS") <(echo "$UNIV_IDS") || true; } | head -10 || true
fi

# ---------------------------------------------------------------------------
# 4. Benchmark + test summary
# ---------------------------------------------------------------------------
bench_report

echo ""
echo "========================================"
echo "  TEST SUMMARY"
echo "========================================"
printf "  Original pipeline : %d circRNAs\n" "${ORIG_CIRCS}"
printf "  Decoupled pipeline: %d circRNAs\n" "${DECOUPLED_CIRCS}"
echo "  PASSED : $PASS"
echo "  FAILED : $FAIL"
echo "========================================"
if [[ $FAIL -eq 0 ]]; then
    echo "  ALL TESTS PASSED"
    exit 0
else
    echo "  SOME TESTS FAILED - see output above"
    exit 1
fi
