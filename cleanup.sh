#!/bin/bash

# Cleanup script for mouse_analysis directory
# Removes old test files, debug outputs, and unused directories

cd /Users/jacobclifford/Desktop/scTE/mouse_analysis

echo "=========================================="
echo "Cleaning up mouse_analysis directory"
echo "=========================================="
echo ""

echo "=== Removing old test/debug scripts ==="
rm -f 05_star_permissive_test.sh
rm -f 06_ultra_basic_star.sh
rm -f 07_star_debug.sh
echo "✓ Removed old test scripts"

echo ""
echo "=== Removing old log and documentation files ==="
rm -f star_debug.log
rm -f test_sequences.fa
rm -f READMEsra.md
rm -f CLAUDE_CODE_INSTRUCTIONS.md
rm -f MEMORY_OPTIMIZATIONS.md
echo "✓ Removed old documentation files (consolidated into README.md)"

echo ""
echo "=== Removing old STAR output from failed attempts ==="
rm -f SRR5068882_no_conda_*.bam
rm -f SRR5068882_no_conda_*.out
rm -rf SRR5068882_no_conda__STARtmp/
echo "✓ Removed old STAR temporary files"

echo ""
echo "=== Removing unused directories ==="
rm -rf docker_workspace/
rm -rf bam/
rm -rf star_output/
rm -rf figures/
echo "✓ Removed unused directories"

echo ""
echo "=== Checking scTE_output directory ==="
if [ -d "scTE_output" ]; then
    SIZE=$(du -sh scTE_output/ | cut -f1)
    echo "scTE_output exists ($SIZE)"
    echo "Keep this? It may contain analysis results."
    echo "To remove manually: rm -rf scTE_output/"
else
    echo "✓ scTE_output does not exist"
fi

echo ""
echo "=== Cleaning up scripts/ folder ==="
echo "Keeping only: run_full_pipeline.sh"

# Remove all scripts except the one we're using
cd scripts/
for file in *.sh; do
    if [ "$file" != "run_full_pipeline.sh" ]; then
        echo "  Removing: $file"
        rm -f "$file"
    fi
done
cd ..

# Remove archive if it exists
if [ -d "scripts/archive" ]; then
    echo "Removing scripts/archive/ directory..."
    rm -rf scripts/archive/
fi

echo "✓ Scripts folder cleaned (only run_full_pipeline.sh remains)"

echo ""
echo "=== Final directory structure ==="
echo ""
echo "KEPT (in git):"
ls -lh Dockerfile docker_*.sh *.md 2>/dev/null
echo ""
ls -lh scripts/*.sh 2>/dev/null

echo ""
echo "DATA DIRECTORIES (kept, in .gitignore):"
du -sh fastq/ reference/ results/ 2>/dev/null

echo ""
echo "=========================================="
echo "Cleanup complete!"
echo "=========================================="
echo ""
echo "To push to GitHub:"
echo "  git push -u origin main"
