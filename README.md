# scTE Mouse Single-Cell RNA-seq Analysis Pipeline

Docker-based pipeline for STAR alignment and transposable element quantification in mouse single-cell RNA-seq data using scTE.

## Overview

This pipeline processes mouse scRNA-seq data (SRR5068882) from C1 Fluidigm platform using STAR alignment and scTE quantification, optimized for 52bp single-end reads. The Docker approach bypasses Apple Silicon and conda compatibility issues encountered with native STAR installations. **scTE is fully integrated and runs automatically** to quantify transposable element expression alongside gene expression.

## Dataset Information

- **SRA ID**: SRR5068882
- **Platform**: C1 Fluidigm single-cell
- **Read Length**: 52bp single-end
- **Total Reads**: ~135 million (6.5GB compressed)
- **Sequencing**: Illumina
- **Species**: Mouse (Mus musculus)
- **Genome**: GRCm38/mm10
- **Annotation**: GENCODE vM21

## Prerequisites

- **Docker Desktop** installed and running
- **~50GB** free disk space
- **16GB RAM** minimum (recommended)
- **Data Files** (if you have them):
  - FASTQ file in `fastq/SRR5068882.fastq.gz`
  - Reference genome in `reference/` (optional - will download if missing)

## Quick Start

### 1. Clone Repository

```bash
git clone https://github.com/jcbclffrd/RNAseqSTAR.git
cd RNAseqSTAR
```

### 2. Place Your Data (Optional)

If you already have the data, place it in the appropriate directories:

```bash
# FASTQ file
mkdir -p fastq
# Copy your SRR5068882.fastq.gz here

# Reference files (optional - will download if not present)
mkdir -p reference
# Copy genome.fa and genes.gtf here if you have them
```

### 3. Build Docker Image

```bash
./docker_build.sh
```

This takes 2-5 minutes and installs STAR and samtools via conda.

### 4. Run Pipeline

```bash
./docker_run_pipeline.sh
```

The pipeline will:
1. Check for/download reference genome and GTF annotation
2. Build STAR genome index with memory-optimized parameters
3. Align reads using STAR (optimized for 52bp reads)
4. Index BAM file with samtools
5. Build scTE index for transposable element annotations
6. Run scTE quantification on aligned reads
7. Display alignment and quantification statistics

**Expected Runtime**: 45-120 minutes for full pipeline (including scTE)

## Pipeline Steps Explained

### Step 1: Data Download (automatic)
- Automatically downloads SRR5068882 using SRA toolkit
- Converts to FASTQ format (~33GB uncompressed)
- Skipped if FASTQ already exists in `data/`

### Step 2: Reference Files (if needed)
- Downloads mouse genome (GRCm38) from Ensembl
- Downloads GENCODE vM21 gene annotations
- **Automatically fixes chromosome naming** (converts "chr1" to "1" format for STAR compatibility)
- Skipped if files already exist in `reference/`

### Step 3: STAR Genome Index
- Builds index with **memory-optimized parameters**:
  - `--genomeSAsparseD 3` - Reduces RAM by using sparser suffix array
  - `--genomeSAindexNbases 12` - Reduces index size for mouse genome
  - `--limitGenomeGenerateRAM 15000000000` - Limits to 15GB during build
  - `--sjdbOverhang 51` - Optimized for 52bp reads (read length - 1)
- Takes 10-30 minutes
- **Peak RAM usage**: ~15GB
- **Disk space**: ~8GB for index

### Step 4: STAR Alignment
- Aligns ~135M reads with parameters optimized for short reads:
  - `--limitBAMsortRAM 4000000000` - 4GB for BAM sorting
  - `--outFilterMultimapNmax 100` - Allows multimapping (important for TEs)
  - `--winAnchorMultimapNmax 100` - Window anchor multimap max
  - `--outFilterMismatchNmax 10` - Max 10 mismatches
  - `--outFilterScoreMinOverLread 0.3` - Min 30% alignment score
  - `--outFilterMatchNmin 15` - Min 15bp matches
- Takes 30-90 minutes depending on system
- **Peak RAM usage**: ~10GB (genome + sorting)
- **Output**: Sorted, coordinate-indexed BAM file

### Step 5: BAM Indexing
- Creates `.bai` index file with samtools
- Shows alignment statistics and file info

### Step 6: scTE Index Building
- Downloads RepeatMasker annotations for mouse (mm10) from UCSC
- Downloads GENCODE gene annotations if needed
- Builds scTE index for efficient TE quantification
- Takes 5-15 minutes
- Skipped if index already exists in `reference/scTE_index/`

### Step 7: scTE Quantification
- Quantifies both gene and transposable element (TE) expression
- Configured for C1 Fluidigm data (no cell barcodes in reads, no UMI)
- Uses multi-threaded processing for speed
- Outputs results in both CSV and HDF5 (.h5ad) format
- Takes 15-30 minutes depending on system

## Memory Optimization Details

This pipeline uses memory-optimized STAR parameters designed for systems with **16GB RAM**:

### Why These Parameters?

Standard STAR genome indexing for mouse requires 30-35GB RAM. Our optimizations reduce this to ~15GB:

| Parameter | Value | Purpose | RAM Saved |
|-----------|-------|---------|-----------|
| `--genomeSAsparseD` | 3 | Sparser suffix array | ~2-3GB |
| `--genomeSAindexNbases` | 12 | Smaller index | ~1-2GB |
| `--limitGenomeGenerateRAM` | 15000000000 | Hard limit | Prevents overflow |
| `--limitBAMsortRAM` | 4000000000 | BAM sort limit | Predictable usage |

### Trade-offs

- **Slightly slower** alignment (negligible for most use cases)
- **Slightly larger** temporary files during BAM sorting
- **No loss** in alignment quality or accuracy

## Output Files

After successful completion, you'll find in `results/`:

```
results/
├── SRR5068882_Aligned.sortedByCoord.out.bam    # Aligned reads (BAM)
├── SRR5068882_Aligned.sortedByCoord.out.bam.bai # BAM index
├── SRR5068882_Log.final.out                     # Alignment statistics
├── SRR5068882_Log.out                           # Detailed log
├── SRR5068882_Log.progress.out                  # Progress during run
└── scTE_output/                                 # scTE quantification results
    ├── scTE_results.csv                         # Gene expression matrix (CSV)
    ├── scTE_results.TE.csv                      # TE expression matrix (CSV)
    └── scTE_results.h5ad                        # Combined data (HDF5/AnnData)
```

### Expected Alignment Statistics

For this C1 Fluidigm dataset (SRR5068882):
- **Input reads**: ~135,000,000
- **Uniquely mapped**: ~86% (excellent for scRNA-seq)
- **Multi-mapped**: ~12%
- **Unmapped**: ~3%

### scTE Output Files

- **scTE_results.csv**: Gene expression counts per cell
- **scTE_results.TE.csv**: Transposable element expression counts per cell
- **scTE_results.h5ad**: Combined gene and TE expression in AnnData format (compatible with Scanpy/Seurat)

## Next Steps: Downstream Analysis

After the pipeline completes, you can analyze the scTE output using standard single-cell analysis tools:

### Using Python (Scanpy)

```python
import scanpy as sc

# Load the scTE results
adata = sc.read_h5ad('results/scTE_output/scTE_results.h5ad')

# The data includes both genes and TEs
print(f"Total features: {adata.n_vars}")
print(f"Total cells: {adata.n_obs}")

# Perform standard scRNA-seq analysis
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable features (genes and TEs)
sc.pp.highly_variable_genes(adata)

# Continue with clustering, UMAP, etc.
```

### Using R (Seurat)

```r
library(Seurat)
library(SeuratDisk)

# Convert h5ad to h5seurat format
Convert("results/scTE_output/scTE_results.h5ad", dest = "h5seurat", overwrite = TRUE)

# Load into Seurat
seurat_obj <- LoadH5Seurat("results/scTE_output/scTE_results.h5seurat")

# Standard Seurat workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
```

## File Structure

```
RNAseqSTAR/
├── Dockerfile                  # Docker configuration
├── docker_build.sh             # Build Docker image
├── docker_run_pipeline.sh      # Run complete pipeline
├── README.md                   # This file
├── .gitignore                  # Git ignore rules
├── scripts/
│   └── run_full_pipeline.sh   # Main pipeline script
├── data/                       # Input FASTQ files (gitignored)
├── reference/                  # Reference genome/GTF (gitignored)
│   └── scTE_index/            # scTE index files (gitignored)
└── results/                    # Output BAM and scTE files (gitignored)
    └── scTE_output/           # scTE quantification results
```

## Troubleshooting

### Docker Issues

**Problem**: `Cannot connect to the Docker daemon`
```bash
# Solution: Start Docker Desktop
# Look for whale icon in menu bar (macOS)
docker info  # Verify it's running
```

**Problem**: Docker build fails
```bash
# Solution: Clean Docker cache
docker system prune -a
# Then rebuild
./docker_build.sh
```

### Memory Issues

**Problem**: System runs out of RAM during genome indexing

**Solutions**:
1. Reduce threads in `scripts/run_full_pipeline.sh`:
   ```bash
   THREADS=4  # Change from 8 to 4
   ```

2. Increase sparsity (further reduce RAM):
   ```bash
   --genomeSAsparseD 4  # Change from 3 to 4
   ```

3. Reduce BAM sort RAM:
   ```bash
   --limitBAMsortRAM 2000000000  # Change from 4GB to 2GB
   ```

### Alignment Issues

**Problem**: 0 input reads (if running natively on macOS)

**Solution**: This is why we use Docker! The Docker approach:
- ✅ Bypasses Apple Silicon compatibility issues
- ✅ Provides clean Linux environment for STAR
- ✅ Avoids conda interference problems
- ✅ Uses pre-compiled conda packages

**Problem**: Very low mapping rate (<30%)

**Possible causes**:
1. Wrong genome version - verify you're using mm10/GRCm38
2. Contamination - check FastQC results
3. Adapter sequences not trimmed

### Disk Space Issues

**Problem**: "No space left on device"

**Solution**: 
```bash
# Check available space
df -h

# Clean up Docker
docker system prune -a

# Remove large intermediate files if needed
rm -rf results/*_STARtmp/
```

## Why Docker?

This pipeline uses Docker to solve several common issues:

| Issue | Native Installation | Docker Solution |
|-------|-------------------|-----------------|
| Apple Silicon compatibility | ❌ STAR fails | ✅ Linux environment |
| Conda interference | ❌ 0 input reads | ✅ Clean environment |
| Dependency conflicts | ❌ Version issues | ✅ Isolated packages |
| Reproducibility | ❌ System-dependent | ✅ Consistent across systems |

## Technical Details

### STAR Parameters Explained

**For 52bp reads (C1 Fluidigm)**:
- `--sjdbOverhang 51` - Optimized for read length (52-1=51)
- Conservative filtering needed for short reads
- Allow multimapping (transposable elements often have multiple copies)

**Memory optimization**:
- Designed for 16GB RAM systems
- Reduces genome index RAM from 30GB to 15GB
- No quality loss, minimal speed impact

### Docker Image

- **Base**: `continuumio/miniconda3:latest`
- **STAR**: v2.7.11b (via conda)
- **samtools**: latest (via conda)
- **Size**: ~2GB after build

## References

- **scTE**: https://github.com/JiekaiLab/scTE
- **STAR aligner**: https://github.com/alexdobin/STAR
- **Dataset**: GSE100058 (He et al., 2017)
- **STAR manual**: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

## Support

For issues specific to this pipeline, please open an issue on GitHub.

For scTE-specific questions, refer to the [scTE documentation](https://github.com/JiekaiLab/scTE).

For STAR-specific questions, refer to the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

## License

This pipeline is provided as-is for research purposes.

---

**Last Updated**: October 2025
