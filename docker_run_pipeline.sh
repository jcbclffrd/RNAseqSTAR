#!/bin/bash

# Run the complete scTE pipeline in Docker container

echo "=========================================="
echo "Running scTE Pipeline in Docker"
echo "=========================================="
echo ""

# Paths on host machine
HOST_BASE="/Users/jacobclifford/Desktop/scTE/mouse_analysis"
HOST_DATA="$HOST_BASE/data"
HOST_REFERENCE="$HOST_BASE/reference"
HOST_RESULTS="$HOST_BASE/results"

# Create directories if they don't exist
mkdir -p "$HOST_DATA" "$HOST_REFERENCE" "$HOST_RESULTS"

echo "Mounting directories:"
echo "  Data: $HOST_DATA"
echo "  Reference: $HOST_REFERENCE"
echo "  Results: $HOST_RESULTS"
echo ""

# Run the Docker container with volume mounts
docker run --rm -it \
    -v "$HOST_DATA:/workspace/data" \
    -v "$HOST_REFERENCE:/workspace/reference" \
    -v "$HOST_RESULTS:/workspace/results" \
    -v "$HOST_BASE/GenomeDir:/workspace/genome_index" \
    --user $(id -u):$(id -g) \
    scte-pipeline:latest \
    /workspace/scripts/run_full_pipeline.sh

echo ""
echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo ""
echo "Results are available in:"
echo "  $HOST_RESULTS"
