#!/usr/bin/env bash
# Build the pyspectrafuse Docker image.
#
# Usage:
#   ./scripts/build_docker.sh                  # local tag + ghcr tag at the current version
#   ./scripts/build_docker.sh 0.0.5            # override version
#   PLATFORM=linux/amd64 ./scripts/build_docker.sh
#   REGISTRY=ghcr.io/bigbio ./scripts/build_docker.sh
#
# Produces:
#   pyspectrafuse:local
#   ${REGISTRY}/pyspectrafuse:${VERSION}
#   ${REGISTRY}/pyspectrafuse:latest
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${HERE}"

VERSION="${1:-$(grep -m1 '^version' pyproject.toml | sed -E 's/.*"([^"]+)".*/\1/')}"
REGISTRY="${REGISTRY:-ghcr.io/bigbio}"
IMAGE="${REGISTRY}/pyspectrafuse"
PLATFORM="${PLATFORM:-linux/amd64}"

echo "Building ${IMAGE}:${VERSION} for ${PLATFORM}"

docker build \
  --platform "${PLATFORM}" \
  -t "pyspectrafuse:local" \
  -t "${IMAGE}:${VERSION}" \
  -t "${IMAGE}:latest" \
  -f Dockerfile \
  .

echo
echo "Tags:"
echo "  pyspectrafuse:local"
echo "  ${IMAGE}:${VERSION}"
echo "  ${IMAGE}:latest"
echo
echo "To push:"
echo "  docker push ${IMAGE}:${VERSION}"
echo "  docker push ${IMAGE}:latest"
