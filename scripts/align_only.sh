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

# Run STAR alignment (output SAM to reduce memory)
echo "Starting STAR alignment (outputting SAM)..."
STAR --runMode alignReads \
     --runThreadN 2 \
     --genomeDir "$GENOME_DIR" \
     --readFilesIn "$INPUT_FASTQ" \
     --readFilesCommand zcat \
     --outFileNamePrefix "$OUTPUT_PREFIX" \
     --outSAMtype SAM \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --sjdbOverhang 51 \
     --outFilterMismatchNmax 10 \
     --outFilterScoreMinOverLread 0.3 \
     --outFilterMatchNmin 15

echo ""
echo "✓ Alignment complete!"

# Convert SAM to sorted BAM
SAM_FILE="${OUTPUT_PREFIX}Aligned.out.sam"
BAM_FILE="${OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam"

if [ -f "$SAM_FILE" ]; then
    echo "Converting SAM to sorted BAM..."
    samtools view -bS "$SAM_FILE" | samtools sort -@ 2 -m 1G -o "$BAM_FILE"

    echo "Indexing BAM file..."
    samtools index "$BAM_FILE"

    # Remove SAM file to save space
    echo "Removing intermediate SAM file..."
    rm "$SAM_FILE"

    echo ""
    echo "=== ALIGNMENT STATISTICS ==="
    cat "${OUTPUT_PREFIX}Log.final.out"
    echo ""

    echo "=== BAM FILE INFO ==="
    ls -lh "$BAM_FILE"
    samtools view -c "$BAM_FILE" | xargs echo "Total alignments:"
else
    echo "ERROR: SAM file not found!"
    exit 1
fi

echo ""
echo "=============================================="
echo "Pipeline Complete!"
echo "=============================================="
