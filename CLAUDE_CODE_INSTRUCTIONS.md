# scTE Pipeline - Docker Setup and Execution Instructions

## Context
We're analyzing transposable elements in mouse single-cell RNA-seq data (SRR5068882) using scTE. We've encountered persistent issues with STAR alignment on Apple Silicon Mac due to conda environment interference. The Docker approach will bypass these macOS-specific issues.

## Current State
- ✅ Data downloaded: `/Users/jacobclifford/Desktop/scTE/mouse_analysis/fastq/SRR5068882.fastq.gz` (6.5GB, 135M reads)
- ✅ Reference genome: `/Users/jacobclifford/Desktop/scTE/mouse_analysis/reference/Mus_musculus.GRCm38.dna.primary_assembly.fa`
- ✅ GTF annotation: `/Users/jacobclifford/Desktop/scTE/mouse_analysis/reference/gencode.vM21.annotation.gtf`
- ✅ STAR genome index: `/Users/jacobclifford/Desktop/scTE/GenomeDir` (built with sjdbOverhang=51 for 52bp reads)
- ✅ Docker Desktop installed and running
- ❌ Native STAR alignment fails with 0 input reads (Apple Silicon + conda issue)

## Data Specifications
- Platform: C1 Fluidigm single-cell
- Read length: 52bp single-end
- Sequencing: Illumina
- STAR parameters needed: sjdbOverhang=51, conservative filtering for short reads

## Task: Build and Run Complete Pipeline in Docker

### Step 1: Create Simple Working Dockerfile

Create `/Users/jacobclifford/Desktop/scTE/mouse_analysis/Dockerfile`:

```dockerfile
# Simple Dockerfile using conda for pre-built binaries
FROM continuumio/miniconda3:latest

# Install STAR and samtools via conda (avoids compilation issues)
RUN conda install -c bioconda -c conda-forge -y \
    star=2.7.11b \
    samtools \
    && conda clean -a -y

# Create workspace directories
RUN mkdir -p /workspace/data \
             /workspace/reference \
             /workspace/genome_index \
             /workspace/results \
             /workspace/scripts

WORKDIR /workspace

# Copy pipeline script
COPY scripts/run_full_pipeline.sh /workspace/scripts/run_full_pipeline.sh
RUN chmod +x /workspace/scripts/run_full_pipeline.sh

# Environment variables
ENV GENOME_DIR=/workspace/genome_index
ENV DATA_DIR=/workspace/data
ENV RESULTS_DIR=/workspace/results
ENV REFERENCE_DIR=/workspace/reference

CMD ["/bin/bash"]
```

### Step 2: Create Pipeline Script

Ensure `/Users/jacobclifford/Desktop/scTE/mouse_analysis/scripts/run_full_pipeline.sh` exists with this content:

```bash
#!/bin/bash
set -e

echo "=============================================="
echo "scTE STAR Alignment in Docker"
echo "=============================================="
echo ""

# Use mounted data (already downloaded)
INPUT_FASTQ="/workspace/data/SRR5068882.fastq.gz"
GENOME_DIR="/workspace/genome_index"
OUTPUT_PREFIX="/workspace/results/SRR5068882_"

# Verify inputs
echo "Checking inputs..."
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: FASTQ not found: $INPUT_FASTQ"
    exit 1
fi

if [ ! -d "$GENOME_DIR" ]; then
    echo "ERROR: Genome index not found: $GENOME_DIR"
    exit 1
fi

echo "✓ Input FASTQ: $INPUT_FASTQ"
echo "✓ Genome index: $GENOME_DIR"
echo "✓ Output prefix: $OUTPUT_PREFIX"
echo ""

# Run STAR alignment
echo "Starting STAR alignment..."
STAR --runMode alignReads \
     --runThreadN 8 \
     --genomeDir "$GENOME_DIR" \
     --readFilesIn "$INPUT_FASTQ" \
     --readFilesCommand zcat \
     --outFileNamePrefix "$OUTPUT_PREFIX" \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 4000000000 \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --sjdbOverhang 51 \
     --outFilterMismatchNmax 10 \
     --outFilterScoreMinOverLread 0.3 \
     --outFilterMatchNmin 15

echo ""
echo "✓ Alignment complete!"

# Index BAM file
BAM_FILE="${OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam"
if [ -f "$BAM_FILE" ]; then
    echo "Indexing BAM file..."
    samtools index "$BAM_FILE"
    
    echo ""
    echo "=== ALIGNMENT STATISTICS ==="
    cat "${OUTPUT_PREFIX}Log.final.out"
    echo ""
    
    echo "=== BAM FILE INFO ==="
    ls -lh "$BAM_FILE"
    samtools view -c "$BAM_FILE" | xargs echo "Total alignments:"
else
    echo "ERROR: BAM file not found!"
    exit 1
fi

echo ""
echo "=============================================="
echo "Pipeline Complete!"
echo "=============================================="
```

### Step 3: Create Docker Build Script

Create `/Users/jacobclifford/Desktop/scTE/mouse_analysis/docker_build.sh`:

```bash
#!/bin/bash

echo "=========================================="
echo "Building scTE Pipeline Docker Image"
echo "=========================================="
echo ""

cd /Users/jacobclifford/Desktop/scTE/mouse_analysis

docker build -t scte-pipeline:latest .

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Docker image built successfully!"
    echo ""
    echo "Next: ./docker_run.sh"
else
    echo ""
    echo "✗ Docker build failed!"
    exit 1
fi
```

### Step 4: Create Docker Run Script

Create `/Users/jacobclifford/Desktop/scTE/mouse_analysis/docker_run.sh`:

```bash
#!/bin/bash

echo "=========================================="
echo "Running STAR Alignment in Docker"
echo "=========================================="
echo ""

# Host paths
HOST_BASE="/Users/jacobclifford/Desktop/scTE/mouse_analysis"
HOST_FASTQ="$HOST_BASE/fastq"
HOST_GENOME="/Users/jacobclifford/Desktop/scTE/GenomeDir"
HOST_RESULTS="$HOST_BASE/results"

# Create results directory
mkdir -p "$HOST_RESULTS"

echo "Mounting directories:"
echo "  FASTQ: $HOST_FASTQ"
echo "  Genome: $HOST_GENOME"
echo "  Results: $HOST_RESULTS"
echo ""

# Run Docker container
docker run --rm -it \
    -v "$HOST_FASTQ:/workspace/data:ro" \
    -v "$HOST_GENOME:/workspace/genome_index:ro" \
    -v "$HOST_RESULTS:/workspace/results" \
    scte-pipeline:latest \
    /workspace/scripts/run_full_pipeline.sh

echo ""
echo "=========================================="
echo "Complete! Results in: $HOST_RESULTS"
echo "=========================================="
```

### Step 5: Execute the Pipeline

Run these commands in order:

```bash
cd /Users/jacobclifford/Desktop/scTE/mouse_analysis

# Make scripts executable
chmod +x docker_build.sh docker_run.sh scripts/run_full_pipeline.sh

# Build Docker image (takes 2-5 minutes)
./docker_build.sh

# Run the pipeline (takes 30-60 minutes)
./docker_run.sh
```

## Expected Output

After successful completion:
- BAM file: `/Users/jacobclifford/Desktop/scTE/mouse_analysis/results/SRR5068882_Aligned.sortedByCoord.out.bam`
- BAM index: `/Users/jacobclifford/Desktop/scTE/mouse_analysis/results/SRR5068882_Aligned.sortedByCoord.out.bam.bai`
- Logs: `/Users/jacobclifford/Desktop/scTE/mouse_analysis/results/SRR5068882_Log.final.out`

## Next Steps After Alignment

Once the BAM file is created, run scTE quantification:

```bash
# Activate scTE environment
cd /Users/jacobclifford/Desktop/scTE
source .venv/bin/activate

# Run scTE (modify parameters based on your cell barcode strategy)
scTE -i mouse_analysis/results/SRR5068882_Aligned.sortedByCoord.out.bam \
     -o mouse_analysis/scTE_output \
     -x mm10 \
     -p 8
```

## Why This Works

1. **Conda in Docker**: Uses pre-compiled STAR binaries (no compilation issues)
2. **Linux environment**: Bypasses Apple Silicon compatibility problems
3. **No conda interference**: Clean environment, no macOS conda conflicts
4. **Volume mounts**: Uses existing data (no re-downloading)
5. **Proven approach**: Uses stable, tested components

## Troubleshooting

If build fails:
- Check Docker Desktop is running: `docker info`
- Check disk space: `df -h`
- Try: `docker system prune -a` to clean up

If run fails:
- Check mounted paths exist
- Verify FASTQ file: `ls -lh /Users/jacobclifford/Desktop/scTE/mouse_analysis/fastq/`
- Verify genome index: `ls -lh /Users/jacobclifford/Desktop/scTE/GenomeDir/`

## Success Criteria

✅ Docker build completes without errors  
✅ STAR reports >0 input reads  
✅ BAM file created with reasonable size (>1GB expected)  
✅ Log shows >50% reads mapped  
✅ Ready for scTE quantification
