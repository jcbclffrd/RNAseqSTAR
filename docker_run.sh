#!/bin/bash

echo "=========================================="
echo "Running STAR Alignment in Docker"
echo "=========================================="
echo ""

# Host paths
HOST_BASE="/Users/jacobclifford/Desktop/scTE/mouse_analysis"
HOST_FASTQ="$HOST_BASE/fastq"
HOST_GENOME="$HOST_BASE/reference/star_genome"
HOST_RESULTS="$HOST_BASE/results"

# Create results directory
mkdir -p "$HOST_RESULTS"

echo "Mounting directories:"
echo "  FASTQ: $HOST_FASTQ"
echo "  Genome: $HOST_GENOME"
echo "  Results: $HOST_RESULTS"
echo ""

# Run Docker container with memory limit
docker run --rm \
    --memory="20g" \
    --memory-swap="20g" \
    -v "$HOST_FASTQ:/workspace/data:ro" \
    -v "$HOST_GENOME:/workspace/genome_index:ro" \
    -v "$HOST_RESULTS:/workspace/results" \
    scte-pipeline:latest \
    /workspace/scripts/align_only.sh

echo ""
echo "=========================================="
echo "Complete! Results in: $HOST_RESULTS"
echo "=========================================="
