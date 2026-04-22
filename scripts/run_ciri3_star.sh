#!/usr/bin/env bash
# =============================================================================
# run_ciri3_star.sh
#
# Production pipeline: CIRI3 four-stage decoupled pipeline for STAR-aligned data.
#
# Edit the "Configuration" section below, then run:
#   bash scripts/run_ciri3_star.sh
#
# Each stage is skipped if its outputs already exist, so the script is safe
# to re-run after a partial failure.
# =============================================================================
set -euo pipefail

# =============================================================================
# Configuration — edit these values for your experiment
# =============================================================================

# STAR output directory suffix for each sample.
# The script expects these directories to exist:
#   STAR_output_${SAMPLE}/Chimeric.out.junction
#   STAR_output_${SAMPLE}/Aligned.out.sam
#   STAR_output_${SAMPLE}/bwa.sam
SAMPLES=(
    "Div_100_S91"
    "Div_101_S92"
    "PARDOS_1_S1"
    "PARDOS_2_S2"
)

# Friendly sample IDs (one per entry in SAMPLES, same order).
# Used as per-sample output folder names and as column headers in the output
# BSJ_Matrix and FSJ_Matrix files.
SAMPLE_IDS=(
    "Control_1"
    "Control_2"
    "Case_1"
    "Case_2"
)

# Directory that contains the STAR_output_* subdirectories.
# Set to ${PWD} if the STAR output dirs are in the current working directory.
STAR_BASE_DIR="${PWD}"

# Reference files
REF_FASTA="/path/to/ref.fa"
GTF_FILE="/path/to/annotation.gtf"

# Number of threads for per-sample stages
THREADS=8

# Path to CIRI3_decoupled.jar
DECOUPLED_JAR="/pastel/Github_scripts/CIRI3/CIRI3_decoupled.jar"

# Output directories for joint stages
UNIVERSE_DIR="${PWD}/universe"
FINALIZE_DIR="${PWD}/finalize"

# =============================================================================
# End of configuration
# =============================================================================

info() { echo "[INFO]  $*"; }
die()  { echo "[ERROR] $*" >&2; exit 1; }

[[ ${#SAMPLES[@]} -eq ${#SAMPLE_IDS[@]} ]] \
    || die "SAMPLES and SAMPLE_IDS must have the same number of entries (got ${#SAMPLES[@]} and ${#SAMPLE_IDS[@]})"

# Locate Java
if [[ -n "${CONDA_PREFIX:-}" ]] && [[ -x "${CONDA_PREFIX}/bin/java" ]]; then
    JAVA_BIN="${CONDA_PREFIX}/bin/java"
else
    JAVA_BIN="$(command -v java)"
fi

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
info "=== Pre-flight checks ==="
[[ -f "${REF_FASTA}"     ]] || die "Reference FASTA not found: ${REF_FASTA}"
[[ -f "${GTF_FILE}"      ]] || die "GTF not found: ${GTF_FILE}"
[[ -x "${JAVA_BIN}"      ]] || die "Java not found: ${JAVA_BIN}"
[[ -f "${DECOUPLED_JAR}" ]] || die "Decoupled jar not found: ${DECOUPLED_JAR}"

for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    STAR_DIR="${STAR_BASE_DIR}/STAR_output_${SAMPLE}"
    [[ -f "${STAR_DIR}/Chimeric.out.junction" ]] || die "Missing: ${STAR_DIR}/Chimeric.out.junction"
    [[ -f "${STAR_DIR}/Aligned.out.sam"       ]] || die "Missing: ${STAR_DIR}/Aligned.out.sam"
    [[ -f "${STAR_DIR}/bwa.sam"               ]] || die "Missing: ${STAR_DIR}/bwa.sam"
done
info "All input files found for ${#SAMPLES[@]} samples."

mkdir -p "${UNIVERSE_DIR}" "${FINALIZE_DIR}"

SCAN1_META_TSV="${UNIVERSE_DIR}/samples_scan1.tsv"
FINALIZE_TSV="${FINALIZE_DIR}/finalize_samples.tsv"

# =============================================================================
# Stage 1: SCAN1 (per sample)
# =============================================================================
info "=== Stage 1: SCAN1 ==="

# Rebuild TSV before the loop so a re-run stays consistent even when
# individual SCAN1 steps are skipped.
> "${SCAN1_META_TSV}"

for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    SAMPLE_ID="${SAMPLE_IDS[$i]}"

    STAR_DIR="${STAR_BASE_DIR}/STAR_output_${SAMPLE}"
    TRIPLE="${STAR_DIR}/Chimeric.out.junction,${STAR_DIR}/Aligned.out.sam,${STAR_DIR}/bwa.sam"
    BWA_SAM="${STAR_DIR}/bwa.sam"

    # Per-sample output folder
    mkdir -p "${PWD}/${SAMPLE_ID}"
    OUT_PREFIX="${PWD}/${SAMPLE_ID}/${SAMPLE_ID}"
    META="${OUT_PREFIX}.scan1_meta"

    if [[ ! -s "${META}" ]]; then
        # Wipe stale BSJ files at both old and new locations before running.
        rm -f "${OUT_PREFIX}BSJ"* "${BWA_SAM}BSJ"*
        info "  SCAN1: ${SAMPLE_ID} (STAR dir: ${SAMPLE})"
        "${JAVA_BIN}" -jar "${DECOUPLED_JAR}" SCAN1 \
            -I  "${TRIPLE}" \
            -O  "${OUT_PREFIX}" \
            -F  "${REF_FASTA}" \
            -A  "${GTF_FILE}" \
            -T  "${THREADS}" \
            -Ma 1 \
            -S  0 \
            2>&1 | tee "${SAMPLE_ID}/scan1.log" \
                 | grep -E "scan|meta|time|Mapped" || true
    else
        info "  [SKIP] SCAN1 already done for ${SAMPLE_ID}"
    fi

    [[ -s "${META}" ]] || die "SCAN1 produced no meta file: ${META}"

    # Always append so TSV is coherent even on a re-run with skipped stages.
    echo -e "${BWA_SAM}\t${META}" >> "${SCAN1_META_TSV}"
done

# =============================================================================
# Stage 2: BUILD_UNIVERSE (joint)
# =============================================================================
info "=== Stage 2: BUILD_UNIVERSE ==="
UNIVERSE_FILE="${UNIVERSE_DIR}/cohort.universe"

if [[ ! -s "${UNIVERSE_FILE}" ]]; then
    "${JAVA_BIN}" -jar "${DECOUPLED_JAR}" BUILD_UNIVERSE \
        -I  "${SCAN1_META_TSV}" \
        -F  "${REF_FASTA}" \
        -O  "${UNIVERSE_DIR}/cohort" \
        2>&1 | tee "${UNIVERSE_DIR}/build_universe.log" \
             | grep -E "Universe|circRNA|time" || true
else
    info "  [SKIP] Universe already exists: ${UNIVERSE_FILE}"
fi

[[ -s "${UNIVERSE_FILE}" ]] || die "BUILD_UNIVERSE produced no universe file: ${UNIVERSE_FILE}"
UNIVERSE_CIRCS=$(grep -c "^chr" "${UNIVERSE_FILE}" || true)
info "Universe: ${UNIVERSE_CIRCS} circRNA candidates."

# =============================================================================
# Stage 3: SCAN2 (per sample)
# =============================================================================
info "=== Stage 3: SCAN2 ==="

# Rebuild TSV before the loop.
> "${FINALIZE_TSV}"

for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    SAMPLE_ID="${SAMPLE_IDS[$i]}"

    STAR_DIR="${STAR_BASE_DIR}/STAR_output_${SAMPLE}"
    TRIPLE="${STAR_DIR}/Chimeric.out.junction,${STAR_DIR}/Aligned.out.sam,${STAR_DIR}/bwa.sam"
    BWA_SAM="${STAR_DIR}/bwa.sam"

    OUT_PREFIX="${PWD}/${SAMPLE_ID}/${SAMPLE_ID}"
    META="${OUT_PREFIX}.scan1_meta"
    FSJ_COUNTS="${OUT_PREFIX}.fsj_counts"

    SPLIT_NUM=$(grep "^fileSplitNum=" "${META}" | cut -d= -f2)

    if [[ ! -s "${FSJ_COUNTS}" ]]; then
        info "  SCAN2: ${SAMPLE_ID}"
        "${JAVA_BIN}" -jar "${DECOUPLED_JAR}" SCAN2 \
            -I  "${TRIPLE}" \
            -CU "${UNIVERSE_FILE}" \
            -SM "${META}" \
            -O  "${OUT_PREFIX}" \
            -F  "${REF_FASTA}" \
            -T  "${THREADS}" \
            -Ma 1 \
            2>&1 | tee "${SAMPLE_ID}/scan2.log" \
                 | grep -E "scan|FSJ|BSJ|time" || true
    else
        info "  [SKIP] SCAN2 already done for ${SAMPLE_ID}"
    fi

    [[ -s "${FSJ_COUNTS}" ]] || die "SCAN2 produced no fsj_counts: ${FSJ_COUNTS}"

    # Col 4 = SAMPLE_ID (column header in output matrices)
    # Col 5 = OUT_PREFIX (BSJ prefix for FINALIZE)
    echo -e "${BWA_SAM}\t${FSJ_COUNTS}\t${SPLIT_NUM}\t${SAMPLE_ID}\t${OUT_PREFIX}" >> "${FINALIZE_TSV}"
done

# =============================================================================
# Stage 4: FINALIZE (joint)
# =============================================================================
info "=== Stage 4: FINALIZE ==="
FINAL_BSJ="${FINALIZE_DIR}/result.BSJ_Matrix"

if [[ ! -s "${FINAL_BSJ}" ]]; then
    "${JAVA_BIN}" -jar "${DECOUPLED_JAR}" FINALIZE \
        -I  "${FINALIZE_TSV}" \
        -CU "${UNIVERSE_FILE}" \
        -F  "${REF_FASTA}" \
        -O  "${FINALIZE_DIR}/result" \
        -A  "${GTF_FILE}" \
        -S  0 \
        2>&1 | tee "${FINALIZE_DIR}/finalize.log" \
             | grep -E "FINALIZE|Summary|Matrix|circRNA|time" || true
else
    info "  [SKIP] FINALIZE already done"
fi

[[ -s "${FINAL_BSJ}" ]]                                  || die "FINALIZE produced no BSJ_Matrix"
[[ -s "${FINALIZE_DIR}/result.FSJ_Matrix" ]]             || die "FINALIZE produced no FSJ_Matrix"

N_CIRCS=$(tail -n +2 "${FINAL_BSJ}" | wc -l)
info "Done. ${N_CIRCS} circRNAs written to:"
info "  ${FINALIZE_DIR}/result.BSJ_Matrix"
info "  ${FINALIZE_DIR}/result.FSJ_Matrix"
