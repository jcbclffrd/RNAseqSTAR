#!/bin/bash

# Master pipeline script for scTE transposable element analysis
# This script runs the complete pipeline from data download to quantification
# Designed to run inside Docker container

set -e  # Exit on any error

echo "=============================================="
echo "scTE Transposable Element Analysis Pipeline"
echo "=============================================="
echo ""

# Configuration
SRR_ID="SRR5068882"
GENOME_URL="ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
GTF_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gtf.gz"
THREADS=8
MEMORY_GB=4

# Directories (using environment variables set in Dockerfile)
DATA_DIR=${DATA_DIR:-/workspace/data}
REFERENCE_DIR=${REFERENCE_DIR:-/workspace/reference}
GENOME_DIR=${GENOME_DIR:-/workspace/genome_index}
RESULTS_DIR=${RESULTS_DIR:-/workspace/results}

echo "Pipeline Configuration:"
echo "  SRR ID: $SRR_ID"
echo "  Threads: $THREADS"
echo "  Memory: ${MEMORY_GB}GB"
echo "  Data directory: $DATA_DIR"
echo "  Reference directory: $REFERENCE_DIR"
echo "  Genome index: $GENOME_DIR"
echo "  Results directory: $RESULTS_DIR"
echo ""

# Create directories
mkdir -p "$DATA_DIR" "$REFERENCE_DIR" "$GENOME_DIR" "$RESULTS_DIR"

# ============================================
# STEP 1: Download SRA data (if not present)
# ============================================
echo "=========================================="
echo "STEP 1: Downloading SRA data"
echo "=========================================="

if [ -f "$DATA_DIR/${SRR_ID}.fastq.gz" ]; then
    echo "✓ FASTQ file already exists: $DATA_DIR/${SRR_ID}.fastq.gz"
elif [ -f "$DATA_DIR/${SRR_ID}.fastq" ]; then
    echo "✓ FASTQ file already exists: $DATA_DIR/${SRR_ID}.fastq"
else
    echo "Downloading $SRR_ID..."
    cd "$DATA_DIR"
    
    # Configure SRA toolkit
    vdb-config --restore-defaults
    
    # Download with prefetch
    prefetch "$SRR_ID" || {
        echo "prefetch failed, trying fastq-dump directly..."
        fastq-dump --gzip --skip-technical --readids --read-filter pass \
                   --dumpbase --split-3 --clip "$SRR_ID"
    }
    
    # Convert to FASTQ if needed
    if [ -d "$SRR_ID" ]; then
        fasterq-dump --outdir "$DATA_DIR" --threads "$THREADS" "$SRR_ID"
        pigz -p "$THREADS" "${SRR_ID}.fastq"
        rm -rf "$SRR_ID"
    fi
    
    echo "✓ Download complete"
fi

echo ""

# ============================================
# STEP 2: Download reference genome and GTF
# ============================================
echo "=========================================="
echo "STEP 2: Downloading reference files"
echo "=========================================="

if [ ! -f "$REFERENCE_DIR/genome.fa" ]; then
    echo "Downloading mouse genome..."
    wget -O "$REFERENCE_DIR/genome.fa.gz" "$GENOME_URL"
    gunzip "$REFERENCE_DIR/genome.fa.gz"
    echo "✓ Genome downloaded"
else
    echo "✓ Genome already exists"
fi

if [ ! -f "$REFERENCE_DIR/genes.gtf" ]; then
    echo "Downloading gene annotations..."
    wget -O "$REFERENCE_DIR/genes.gtf.gz" "$GTF_URL"
    gunzip "$REFERENCE_DIR/genes.gtf.gz"
    echo "✓ GTF downloaded"
else
    echo "✓ GTF already exists"
fi

echo ""

# ============================================
# STEP 3: Build STAR genome index
# ============================================
echo "=========================================="
echo "STEP 3: Building STAR genome index"
echo "=========================================="

if [ -f "$GENOME_DIR/SA" ]; then
    echo "✓ STAR genome index already exists"
else
    echo "Building STAR index with memory-optimized parameters..."
    echo "  --genomeSAsparseD 3 (reduces RAM usage)"
    echo "  --genomeSAindexNbases 12 (reduces index size)"
    echo "  --limitGenomeGenerateRAM 15000000000 (15GB limit)"
    echo "This takes 10-30 minutes..."
    
    STAR --runMode genomeGenerate \
         --runThreadN "$THREADS" \
         --genomeDir "$GENOME_DIR" \
         --genomeFastaFiles "$REFERENCE_DIR/genome.fa" \
         --sjdbGTFfile "$REFERENCE_DIR/genes.gtf" \
         --sjdbOverhang 51 \
         --genomeSAsparseD 3 \
         --genomeSAindexNbases 12 \
         --limitGenomeGenerateRAM 15000000000
    
    echo "✓ Genome index built"
fi

echo ""

# ============================================
# STEP 4: STAR alignment
# ============================================
echo "=========================================="
echo "STEP 4: STAR alignment"
echo "=========================================="

# Determine input file
if [ -f "$DATA_DIR/${SRR_ID}.fastq.gz" ]; then
    INPUT_FILE="$DATA_DIR/${SRR_ID}.fastq.gz"
    READ_FILES_COMMAND="zcat"
elif [ -f "$DATA_DIR/${SRR_ID}.fastq" ]; then
    INPUT_FILE="$DATA_DIR/${SRR_ID}.fastq"
    READ_FILES_COMMAND="cat"
else
    echo "ERROR: FASTQ file not found!"
    exit 1
fi

OUTPUT_PREFIX="$RESULTS_DIR/${SRR_ID}_"

echo "Input: $INPUT_FILE"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""
echo "Running STAR with memory-optimized parameters for 52bp reads..."
echo "  --limitBAMsortRAM 4000000000 (4GB for BAM sorting)"
echo ""

STAR --runMode alignReads \
     --runThreadN "$THREADS" \
     --genomeDir "$GENOME_DIR" \
     --readFilesIn "$INPUT_FILE" \
     --readFilesCommand "$READ_FILES_COMMAND" \
     --outFileNamePrefix "$OUTPUT_PREFIX" \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 4000000000 \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --sjdbOverhang 51 \
     --outFilterMismatchNmax 10 \
     --outFilterScoreMinOverLread 0.3 \
     --outFilterMatchNmin 15

echo "✓ Alignment complete"
echo ""

# ============================================
# STEP 5: Index BAM file
# ============================================
echo "=========================================="
echo "STEP 5: Indexing BAM file"
echo "=========================================="

BAM_FILE="${OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam"

if [ -f "$BAM_FILE" ]; then
    echo "Indexing BAM file..."
    samtools index "$BAM_FILE"
    echo "✓ BAM indexed"
    
    # Show alignment statistics
    echo ""
    echo "=== ALIGNMENT STATISTICS ==="
    cat "${OUTPUT_PREFIX}Log.final.out"
    echo ""
    
    # Show BAM stats
    echo "=== BAM FILE INFO ==="
    ls -lh "$BAM_FILE"
    echo "Total alignments: $(samtools view -c "$BAM_FILE")"
    echo ""
else
    echo "ERROR: BAM file not found: $BAM_FILE"
    exit 1
fi

# ============================================
# PIPELINE COMPLETE
# ============================================
echo "=============================================="
echo "Pipeline Complete!"
echo "=============================================="
echo ""
echo "Output files:"
echo "  BAM: $BAM_FILE"
echo "  BAI: ${BAM_FILE}.bai"
echo "  Logs: ${OUTPUT_PREFIX}Log.*.out"
echo ""
echo "Next step: Run scTE quantification"
echo "  scTE -i $BAM_FILE -o $RESULTS_DIR/scTE_output -x mm10 -CB <cell_barcode_tag>"
echo ""
