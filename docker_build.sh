#!/bin/bash

# Build the Docker image for scTE pipeline

echo "=========================================="
echo "Building scTE Pipeline Docker Image"
echo "=========================================="
echo ""

cd /Users/jacobclifford/Desktop/scTE/mouse_analysis

# Build the Docker image
docker build -t scte-pipeline:latest .

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Docker image built successfully!"
    echo ""
    echo "Image: scte-pipeline:latest"
    echo ""
    echo "Next step: Run the pipeline with:"
    echo "  ./docker_run_pipeline.sh"
else
    echo ""
    echo "✗ Docker build failed!"
    exit 1
fi
