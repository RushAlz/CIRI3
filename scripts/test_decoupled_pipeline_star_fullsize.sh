#!/usr/bin/env bash
# =============================================================================
# test_decoupled_pipeline_star_fullsize.sh
#
# Runs the original joint (-W 1, -Ma 1) STAR CIRI3 pipeline AND the four-stage
# decoupled pipeline (SCAN1 -> BUILD_UNIVERSE -> SCAN2 -> FINALIZE) on
# full-size STAR data, then verifies that both pipelines produce identical
# BSJ and FSJ matrices.
#
# The original (joint) pipeline runs the stock CIRI3_Java_18.0.1.jar so the
# comparison is against the published ground truth. The decoupled pipeline
# runs CIRI3_decoupled.jar, the jar built from this repo's src/ tree.
#
# Assumes alignment is already complete (STAR + BWA-MEM re-alignment of
# unmapped reads) and STAR output directories exist under DATA_DIR.
#
# Each stage is skipped if its outputs already exist, so the script is safe
# to re-run after a partial failure.
#
# Usage:
#   bash scripts/test_decoupled_pipeline_star_fullsize.sh [options]
#
#   --threads N          threads for each per-sample stage (default: 8)
#   --output-dir D       output directory (default: DATA_DIR/decoupled_comparison)
#   --intron             run with intron mode (-It 1) in both pipelines
#   --use-current-joint  run the joint pipeline from CIRI3_decoupled.jar (this
#                        repo's current source) instead of CIRI3_Java_18.0.1.jar.
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

# Friendly sample IDs used as per-sample output folder names and as column
# headers in BSJ_Matrix / FSJ_Matrix. Must be in the same order as SAMPLES.
SAMPLE_IDS=(
    "Div_100_S91"
    "Div_101_S92"
    "PARDOS_1_S1"
    "PARDOS_2_S2"
)

# ---------------------------------------------------------------------------
# Derived paths
# ---------------------------------------------------------------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

ORIGINAL_JAR="${REPO_ROOT}/CIRI3_Java_18.0.1.jar"       # published ground truth
DECOUPLED_JAR="${REPO_ROOT}/CIRI3_decoupled.jar"       # built from this repo

[[ -z "$OUT_ROOT" ]] && OUT_ROOT="${DATA_DIR}/decoupled_comparison"

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
    JOINT_LABEL="CIRI3_Java_18.0.1.jar (published, -W 1)"
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
    STAR_DIR="${DATA_DIR}/STAR_output_${S}"
    [[ -f "${STAR_DIR}/Chimeric.out.junction" ]] || die "Missing: ${STAR_DIR}/Chimeric.out.junction"
    [[ -f "${STAR_DIR}/Aligned.out.sam"       ]] || die "Missing: ${STAR_DIR}/Aligned.out.sam"
    [[ -f "${STAR_DIR}/bwa.sam"               ]] || die "Missing: ${STAR_DIR}/bwa.sam"
done
info "All input files found for ${#SAMPLES[@]} samples."
info "Output directory: ${OUT_ROOT}"
info "Joint pipeline: ${JOINT_LABEL}"
info "Decoupled jar:  ${DECOUPLED_JAR}"

mkdir -p "$ORIG_DIR" "$SCAN1_DIR" "$UNIVERSE_DIR" "$SCAN2_DIR" "$FINALIZE_DIR" "$BENCH_DIR"

# ---------------------------------------------------------------------------
# 1. ORIGINAL pipeline (-W 1, STAR) from CIRI3_Java_18.0.1.jar
# ---------------------------------------------------------------------------
info "=== Stage 0: ORIGINAL pipeline (CIRI3_Java_18.0.1.jar, -W 1, -Ma 1) ==="

ORIG_BSJ="${ORIG_DIR}/result.BSJ_Matrix"
if [[ -s "$ORIG_BSJ" ]]; then
    info "  [SKIP] Original results exist at ${ORIG_DIR}/result"
else
    ORIG_TSV="${ORIG_DIR}/samples.tsv"
    > "$ORIG_TSV"
    for S in "${SAMPLES[@]}"; do
        STAR_DIR="${DATA_DIR}/STAR_output_${S}"
        # Wipe any BSJ files left over from a previous test run; joint opens
        # BSJ{threadNum} in APPEND mode, so stale records would accumulate.
        rm -f "${STAR_DIR}/bwa.sam"BSJ*
        echo "${STAR_DIR}/Chimeric.out.junction,${STAR_DIR}/Aligned.out.sam,${STAR_DIR}/bwa.sam" \
            >> "$ORIG_TSV"
    done

    bench_run "00_original_joint" \
        "${JAVA_ORIG[@]}" \
            -I "${ORIG_TSV}" \
            -O "${ORIG_DIR}/result" \
            -F "${REF_FA}" \
            -A "${GTF_FILE}" \
            -W 1 -Ma 1 -T "${THREADS}" -S 2 "${INTRON_FLAG[@]}" \
        2>&1 | tee "${ORIG_DIR}/run.log" \
        | grep -E "CIRI3|scan|completed|circRNA|Mapped|time|Exception|Error|^\t?at |DIAG" || true
fi

if [[ ! -s "${ORIG_DIR}/result.BSJ_Matrix" ]]; then
    echo "[ERROR] Original pipeline produced no BSJ_Matrix — aborting." >&2
    echo "        Inspect: ${ORIG_DIR}/run.log" >&2
    tail -40 "${ORIG_DIR}/run.log" >&2 || true
    exit 1
fi

# Snapshot joint BSJ files for later direct diff against decoupled BSJ files.
# Only effective when the user runs with CIRI3_KEEP_BSJ=1 (which tells current-
# source MutFileSTARTest to skip its internal BSJ deletion). If the env-var
# wasn't set, the joint BSJ files are already gone and we just skip.
JOINT_BSJ_SNAPSHOT="${OUT_ROOT}/joint_bsj_snapshot"
if [[ -n "${CIRI3_KEEP_BSJ:-}" ]]; then
    mkdir -p "${JOINT_BSJ_SNAPSHOT}"
    for S in "${SAMPLES[@]}"; do
        STAR_DIR="${DATA_DIR}/STAR_output_${S}"
        BWA_SAM="${STAR_DIR}/bwa.sam"
        mkdir -p "${JOINT_BSJ_SNAPSHOT}/${S}"
        # Move (not copy) so they don't interfere with decoupled SCAN1's pre-clean.
        for f in "${BWA_SAM}BSJ"*; do
            [[ -f "$f" ]] || continue
            mv "$f" "${JOINT_BSJ_SNAPSHOT}/${S}/"
        done
    done
    info "Joint BSJ files snapshotted to ${JOINT_BSJ_SNAPSHOT}/"
fi

check_exists "Original BSJ_Matrix" "${ORIG_DIR}/result.BSJ_Matrix"
check_exists "Original FSJ_Matrix" "${ORIG_DIR}/result.FSJ_Matrix"
ORIG_CIRCS=$(tail -n +2 "${ORIG_DIR}/result.BSJ_Matrix" | wc -l)
info "Original pipeline: ${ORIG_CIRCS} circRNAs."

# ---------------------------------------------------------------------------
# 2. DECOUPLED pipeline (CIRI3_decoupled.jar, STAR)
# ---------------------------------------------------------------------------
info "=== Decoupled pipeline (CIRI3_decoupled.jar, STAR) ==="

> "$SCAN1_META_TSV"
> "$FINALIZE_TSV"

# --- Stage 1: SCAN1 (per sample) ---
info "--- Stage 1: SCAN1 ---"
SCAN1_IDX=0
for i in "${!SAMPLES[@]}"; do
    S="${SAMPLES[$i]}"
    SAMPLE_ID="${SAMPLE_IDS[$i]}"
    SCAN1_IDX=$((SCAN1_IDX+1))
    STAR_DIR="${DATA_DIR}/STAR_output_${S}"
    TRIPLE="${STAR_DIR}/Chimeric.out.junction,${STAR_DIR}/Aligned.out.sam,${STAR_DIR}/bwa.sam"
    BWA_SAM="${STAR_DIR}/bwa.sam"
    mkdir -p "${SCAN1_DIR}/${SAMPLE_ID}"
    OUT_PREFIX="${SCAN1_DIR}/${SAMPLE_ID}/${SAMPLE_ID}"
    META="${OUT_PREFIX}.scan1_meta"

    if [[ -s "$META" ]]; then
        info "  [SKIP] SCAN1 already done for ${SAMPLE_ID}"
    else
        # Wipe stale BSJ files at both new per-sample location and old
        # bwa.sam-adjacent location so SCAN2's BSJ enumeration is clean.
        rm -f "${OUT_PREFIX}BSJ"* "${BWA_SAM}BSJ"*
        info "  SCAN1: ${SAMPLE_ID} (STAR dir: ${S})"
        scan1_stdout="${BENCH_DIR}/scan1_${SAMPLE_ID}.stdout"
        bench_run "$(printf '10_scan1_%02d_%s' "${SCAN1_IDX}" "${SAMPLE_ID}")" \
            "${JAVA_NEW[@]}" SCAN1 \
                -I "${SCAN1_PAIR}" \
                -O "${OUT_PREFIX}" \
                -F "${REF_FA}" \
                -A "${GTF_FILE}" \
                -T "${THREADS}" -Ma 1 -S 2 "${INTRON_FLAG[@]}" \
            2>&1 | tee "${scan1_stdout}" | grep -iE "SCAN1|scan.*completed|meta|time|Mapped|Exception|Error|BrokenBarrier" || true
        # If meta still missing, surface the tail of the captured output and the Java log
        if [[ ! -s "$META" ]]; then
            echo "[DEBUG] Last 30 lines of SCAN1 stdout (${scan1_stdout}):"
            tail -30 "${scan1_stdout}" 2>/dev/null || echo "  (no stdout captured)"
            java_log="${OUT_PREFIX}.log"
            if [[ -s "$java_log" ]]; then
                echo "[DEBUG] Last 30 lines of Java log (${java_log}):"
                tail -30 "${java_log}"
            else
                echo "[DEBUG] Java log not found: ${java_log}"
            fi
        fi
    fi

    check_exists "SCAN1 meta (${SAMPLE_ID})" "${META}"

    # Guard: if meta is missing (SCAN1 failed), skip BSJ checks for this sample.
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

    # 1-column meta list consumed by MAKE_BSJ_LIST in the BUILD_UNIVERSE step below
    echo "$META" >> "$SCAN1_META_TSV"
done

# --- Stage 2: BUILD_UNIVERSE (via MAKE_BSJ_LIST + -IB) ---
info "--- Stage 2: BUILD_UNIVERSE ---"
UNIVERSE_FILE="${UNIVERSE_DIR}/cohort.universe"
BSJ_LIST_TSV="${UNIVERSE_DIR}/cohort.bsj_list.tsv"
if [[ -s "$UNIVERSE_FILE" ]]; then
    info "  [SKIP] Universe already exists"
else
    # Convert scan1_meta files to the direct BSJ-list format, then build universe.
    # MAKE_BSJ_LIST reads bsjPrefix/fileSplitNum/readLen from each meta file;
    # no stale VM-local paths are propagated into the universe build step.
    bench_run "19_make_bsj_list" \
        "${JAVA_NEW[@]}" MAKE_BSJ_LIST \
            -I "${SCAN1_META_TSV}" \
            -O "${UNIVERSE_DIR}/cohort" \
        2>&1 | grep -E "MAKE_BSJ|rows|WARNING" || true
    check_exists "BSJ list for BUILD_UNIVERSE" "${BSJ_LIST_TSV}"

    bench_run "20_build_universe" \
        "${JAVA_NEW[@]}" BUILD_UNIVERSE \
            -IB "${BSJ_LIST_TSV}" \
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
for i in "${!SAMPLES[@]}"; do
    S="${SAMPLES[$i]}"
    SAMPLE_ID="${SAMPLE_IDS[$i]}"
    SCAN2_IDX=$((SCAN2_IDX+1))
    STAR_DIR="${DATA_DIR}/STAR_output_${S}"
    TRIPLE="${STAR_DIR}/Chimeric.out.junction,${STAR_DIR}/Aligned.out.sam,${STAR_DIR}/bwa.sam"
    BWA_SAM="${STAR_DIR}/bwa.sam"
    SCAN1_OUT_PREFIX="${SCAN1_DIR}/${SAMPLE_ID}/${SAMPLE_ID}"
    META="${SCAN1_OUT_PREFIX}.scan1_meta"
    SPLIT_NUM=$(grep "^fileSplitNum=" "${META}" | cut -d= -f2)
    mkdir -p "${SCAN2_DIR}/${SAMPLE_ID}"
    OUT_PREFIX="${SCAN2_DIR}/${SAMPLE_ID}/${SAMPLE_ID}"
    FSJ_COUNTS="${OUT_PREFIX}.fsj_counts"

    if [[ -s "$FSJ_COUNTS" ]]; then
        info "  [SKIP] SCAN2 already done for ${SAMPLE_ID}"
    else
        info "  SCAN2: ${SAMPLE_ID}"
        bench_run "$(printf '30_scan2_%02d_%s' "${SCAN2_IDX}" "${SAMPLE_ID}")" \
            "${JAVA_NEW[@]}" SCAN2 \
                -I "${TRIPLE}" \
                -CU "${UNIVERSE_FILE}" \
                -SM "${META}" \
                -O "${OUT_PREFIX}" \
                -F "${REF_FA}" \
                -T "${THREADS}" -Ma 1 "${INTRON_FLAG[@]}" \
            2>&1 | grep -E "scan|FSJ|BSJ|time" || true
    fi
    check_exists "SCAN2 FSJ counts (${SAMPLE_ID})" "${FSJ_COUNTS}"

    # 4-column direct BSJ list format for FINALIZE -IB:
    #   bsjPrefix <TAB> fileSplitNum <TAB> fsjCountsFile <TAB> sampleName
    echo -e "${SCAN1_OUT_PREFIX}\t${SPLIT_NUM}\t${FSJ_COUNTS}\t${SAMPLE_ID}" >> "$FINALIZE_TSV"
done

# --- Stage 4: FINALIZE ---
info "--- Stage 4: FINALIZE ---"
FINAL_BSJ="${FINALIZE_DIR}/result.BSJ_Matrix"
if [[ -s "$FINAL_BSJ" ]]; then
    info "  [SKIP] FINALIZE already done"
else
    bench_run "40_finalize" \
        "${JAVA_NEW[@]}" FINALIZE \
            -IB "${FINALIZE_TSV}" \
            -CU "${UNIVERSE_FILE}" \
            -F "${REF_FA}" \
            -O "${FINALIZE_DIR}/result" \
            -A "${GTF_FILE}" \
            -S 2 "${INTRON_FLAG[@]}" \
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

# Per-sample, per-BSJ-file diff between joint's snapshotted BSJ files and the
# decoupled run's BSJ files. Only active when CIRI3_KEEP_BSJ was set for joint.
if [[ -d "${JOINT_BSJ_SNAPSHOT:-}" ]]; then
    info "=== Per-BSJ-file joint vs decoupled diff ==="
    for idx in "${!SAMPLES[@]}"; do
        S="${SAMPLES[$idx]}"
        SAMPLE_ID="${SAMPLE_IDS[$idx]}"
        STAR_DIR="${DATA_DIR}/STAR_output_${S}"
        BWA_SAM="${STAR_DIR}/bwa.sam"
        DECOUPLED_BSJ_PREFIX="${SCAN1_DIR}/${SAMPLE_ID}/${SAMPLE_ID}"
        for i in 1 2 3 4 5; do
            joint_f="${JOINT_BSJ_SNAPSHOT}/${S}/$(basename "$BWA_SAM")BSJ${i}"
            dec_f="${DECOUPLED_BSJ_PREFIX}BSJ${i}"
            [[ -f "$joint_f" && -f "$dec_f" ]] || continue
            joint_sz=$(wc -l < "$joint_f")
            dec_sz=$(wc -l < "$dec_f")
            if [[ "$joint_sz" == "$dec_sz" ]] && cmp -s \
                   <(sort "$joint_f") <(sort "$dec_f"); then
                info "  ${S} BSJ${i}: IDENTICAL (${joint_sz} lines)"
            else
                fail "  ${S} BSJ${i}: joint=${joint_sz} lines, decoupled=${dec_sz} lines"
                if [[ "$i" == "2" ]]; then
                    echo "    lines in joint only (first 5):"
                    { comm -23 <(sort "$joint_f") <(sort "$dec_f") || true; } \
                        | head -5 | sed 's/^/      /' || true
                    echo "    lines in decoupled only (first 5):"
                    { comm -13 <(sort "$joint_f") <(sort "$dec_f") || true; } \
                        | head -5 | sed 's/^/      /' || true
                fi
            fi
        done
    done
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
