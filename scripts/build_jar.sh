#!/usr/bin/env bash
# =============================================================================
# build_jar.sh
#
# Compiles the CIRI3 sources and packages them (plus htsjdk) into
# CIRI3_decoupled.jar at the repository root. This is the jar that the
# decoupled pipeline scripts run.
#
# Usage:
#   bash scripts/build_jar.sh
# =============================================================================
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="${REPO_ROOT}/bin"
LIB_DIR="${REPO_ROOT}/lib"
SRC_DIR="${REPO_ROOT}/src"
OUT_JAR="${REPO_ROOT}/CIRI3_decoupled.jar"

HTSJDK="${LIB_DIR}/htsjdk-3.0.4.jar"
CLASSPATH="${BIN_DIR}:${HTSJDK}"

if [[ -n "${CONDA_PREFIX:-}" ]] && [[ -x "${CONDA_PREFIX}/bin/javac" ]]; then
    JAVAC_BIN="${CONDA_PREFIX}/bin/javac"
    JAR_BIN="${CONDA_PREFIX}/bin/jar"
else
    JAVAC_BIN="$(command -v javac)"
    JAR_BIN="$(command -v jar)"
fi

[[ -x "$JAVAC_BIN" ]] || { echo "javac not found"; exit 1; }
[[ -x "$JAR_BIN"   ]] || { echo "jar not found"  ; exit 1; }

echo "[INFO] Compiling sources to ${BIN_DIR} (target: Java 8)..."
mkdir -p "${BIN_DIR}"
find "${SRC_DIR}" -name "*.java" > /tmp/ciri3_sources_$$.txt
"${JAVAC_BIN}" -source 8 -target 8 -cp "${CLASSPATH}" -d "${BIN_DIR}" \
    @/tmp/ciri3_sources_$$.txt 2>&1 \
    | grep -v "^\(warning\|Note\)" || true
rm -f /tmp/ciri3_sources_$$.txt

echo "[INFO] Staging jar contents..."
STAGE_DIR="$(mktemp -d)"
trap 'rm -rf "${STAGE_DIR}"' EXIT
cp -r "${BIN_DIR}/." "${STAGE_DIR}/"
( cd "${STAGE_DIR}" && unzip -oq "${HTSJDK}" -x 'META-INF/*' )
mkdir -p "${STAGE_DIR}/META-INF"
cat > "${STAGE_DIR}/META-INF/MANIFEST.MF" <<'EOF'
Manifest-Version: 1.0
Class-Path: .
Main-Class: com.zx.test.TestParameters
EOF

echo "[INFO] Packaging ${OUT_JAR}..."
( cd "${STAGE_DIR}" && "${JAR_BIN}" cfm "${OUT_JAR}" META-INF/MANIFEST.MF . )

echo "[INFO] Done. Built: ${OUT_JAR}"
ls -la "${OUT_JAR}"
