#!/usr/bin/env bash
# Build a Singularity SIF from the pyspectrafuse Docker image.
#
# Two modes:
#   local  — build the Docker image locally, then convert via docker-daemon:
#            (needs Docker + singularity/apptainer on the same host)
#   remote — pull the published image from ghcr.io into a SIF directly
#            (needs only singularity/apptainer; works on HPC login nodes)
#
# Usage:
#   ./scripts/build_singularity.sh                         # local, current version
#   ./scripts/build_singularity.sh 0.0.5                   # local, override version
#   MODE=remote ./scripts/build_singularity.sh             # pull from ghcr.io
#   MODE=remote OUT_DIR=/hps/nobackup/.../singularity \
#       ./scripts/build_singularity.sh                     # into EBI Codon cache
#
# Output SIF filename follows the nf-core/bigbio convention used by the
# EBI Codon profile: ghcr.io-bigbio-pyspectrafuse-<version>.sif
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${HERE}"

VERSION="${1:-$(grep -m1 '^version' pyproject.toml | sed -E 's/.*"([^"]+)".*/\1/')}"
REGISTRY="${REGISTRY:-ghcr.io/bigbio}"
IMAGE="${REGISTRY}/pyspectrafuse"
MODE="${MODE:-local}"
OUT_DIR="${OUT_DIR:-.}"
SIF_NAME="$(echo "${IMAGE}:${VERSION}" | tr ':/' '--').sif"   # e.g. ghcr.io-bigbio-pyspectrafuse-0.0.4.sif
SIF_PATH="${OUT_DIR}/${SIF_NAME}"

mkdir -p "${OUT_DIR}"

# Pick singularity / apptainer
SINGULARITY="${SINGULARITY:-}"
if [[ -z "${SINGULARITY}" ]]; then
  if command -v apptainer >/dev/null 2>&1; then
    SINGULARITY=apptainer
  elif command -v singularity >/dev/null 2>&1; then
    SINGULARITY=singularity
  else
    echo "ERROR: neither apptainer nor singularity found on PATH." >&2
    exit 1
  fi
fi
echo "Using ${SINGULARITY}"

case "${MODE}" in
  local)
    echo "Mode: local — building Docker image then converting to SIF"
    "${HERE}/scripts/build_docker.sh" "${VERSION}"
    echo "Converting docker-daemon://${IMAGE}:${VERSION} → ${SIF_PATH}"
    "${SINGULARITY}" build --force "${SIF_PATH}" "docker-daemon://${IMAGE}:${VERSION}"
    ;;
  remote)
    echo "Mode: remote — pulling ${IMAGE}:${VERSION} directly into SIF"
    "${SINGULARITY}" build --force "${SIF_PATH}" "docker://${IMAGE}:${VERSION}"
    ;;
  *)
    echo "ERROR: MODE must be 'local' or 'remote' (got: ${MODE})" >&2
    exit 1
    ;;
esac

echo
echo "SIF: ${SIF_PATH}"
"${SINGULARITY}" inspect "${SIF_PATH}" || true
