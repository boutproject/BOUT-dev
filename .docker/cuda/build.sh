#!/bin/bash
# Build script for BOUT++ CUDA Docker container

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_ROOT="$( cd "$SCRIPT_DIR/../.." && pwd )"

# Default values
IMAGE_NAME="bout-cuda"
IMAGE_TAG="latest"
BOUT_URL="https://github.com/boutproject/BOUT-dev.git"
BOUT_COMMIT="main"
CUDA_ARCH="80"  # Default to Ampere (A100)
BUILD_ARGS=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -n|--name)
      IMAGE_NAME="$2"
      shift 2
      ;;
    -t|--tag)
      IMAGE_TAG="$2"
      shift 2
      ;;
    --url)
      BOUT_URL="$2"
      shift 2
      ;;
    --commit)
      BOUT_COMMIT="$2"
      shift 2
      ;;
    --cuda-arch)
      CUDA_ARCH="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 [OPTIONS]"
      echo ""
      echo "Options:"
      echo "  -n, --name NAME        Image name (default: bout-cuda)"
      echo "  -t, --tag TAG          Image tag (default: latest)"
      echo "  --url URL              BOUT++ repository URL"
      echo "  --commit COMMIT        BOUT++ commit/branch to build (default: main)"
      echo "  --cuda-arch ARCH       CUDA architecture (default: 80 for Ampere)"
      echo "                         70=Volta, 75=Turing, 80=Ampere, 86=Ampere(mobile)"
      echo "                         89=Ada, 90=Hopper"
      echo "  -h, --help             Show this help message"
      echo ""
      echo "Examples:"
      echo "  $0                                    # Build with defaults"
      echo "  $0 --cuda-arch 70                     # Build for Volta (V100)"
      echo "  $0 --commit v5.2.0                    # Build specific version"
      echo "  $0 --name my-bout --tag dev           # Custom image name and tag"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use -h or --help for usage information"
      exit 1
      ;;
  esac
done

# Update Dockerfile with CUDA architecture if not default
if [ "$CUDA_ARCH" != "80" ]; then
  echo "Note: Building for CUDA architecture $CUDA_ARCH"
  echo "      You may need to edit the Dockerfile to change cuda_arch in Spack commands"
fi

BUILD_ARGS="--build-arg BOUT_URL=${BOUT_URL}"
BUILD_ARGS="${BUILD_ARGS} --build-arg BOUT_COMMIT=${BOUT_COMMIT}"

echo "============================================"
echo "Building BOUT++ CUDA Docker Image"
echo "============================================"
echo "Image name:    ${IMAGE_NAME}:${IMAGE_TAG}"
echo "BOUT++ URL:    ${BOUT_URL}"
echo "BOUT++ commit: ${BOUT_COMMIT}"
echo "CUDA arch:     ${CUDA_ARCH}"
echo "Build context: ${SCRIPT_DIR}"
echo "============================================"
echo ""

# Build the image
cd "$SCRIPT_DIR"

docker build \
  ${BUILD_ARGS} \
  -t "${IMAGE_NAME}:${IMAGE_TAG}" \
  -f Dockerfile \
  .

echo ""
echo "============================================"
echo "Build complete!"
echo "============================================"
echo "Run with: docker run --gpus all -it ${IMAGE_NAME}:${IMAGE_TAG}"
echo ""
echo "To test CUDA:"
echo "  docker run --gpus all -it ${IMAGE_NAME}:${IMAGE_TAG} bash -c 'nvidia-smi'"
echo ""
