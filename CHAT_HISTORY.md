# RNAseqSTAR Project Chat History

## Project Overview
Development of a Docker-based pipeline for STAR alignment and scTE transposable element quantification in mouse single-cell RNA-seq data.

## Key Timeline & Decisions

### Initial Problem
- **Issue**: STAR alignment failing on Apple Silicon Mac with 0 input reads
- **Root Cause**: Conda environment interference + Apple Silicon compatibility issues
- **Dataset**: SRR5068882 (C1 Fluidigm, 52bp single-end, ~135M reads)

### Technical Challenges Solved

#### 1. STAR Compatibility Issues
- **Problem**: Native STAR installation failing with 0 reads despite valid FASTQ
- **Attempts**: Multiple homebrew/conda installations, parameter adjustments
- **Discovery**: Conda base environment was interfering with STAR on Apple Silicon
- **Solution**: Docker approach to bypass macOS-specific issues

#### 2. Memory Optimization
- **Challenge**: Standard STAR genome indexing requires 30-35GB RAM
- **Solution**: Memory-optimized parameters for 16GB systems:
  - `--genomeSAsparseD 3` (reduces RAM by using sparser suffix array)
  - `--genomeSAindexNbases 12` (reduces index size)
  - `--limitGenomeGenerateRAM 15000000000` (15GB limit)
  - `--limitBAMsortRAM 4000000000` (4GB for BAM sorting)

#### 3. Read Length Optimization
- **Specification**: 52bp single-end reads from C1 Fluidigm
- **Parameters**: `--sjdbOverhang 51` (read length - 1)
- **Filtering**: Conservative parameters for short reads

### Architecture Evolution

#### Phase 1: Native Attempts
- Multiple STAR installation methods
- Conda environment troubleshooting
- Parameter optimization attempts
- **Result**: Persistent 0 input reads issue

#### Phase 2: Docker Transition
- Created multiple Docker scripts (03, 04, 05, 06 series)
- Dockerfile compilation issues with building STAR from source
- **Result**: Overly complex multi-script approach

#### Phase 3: Simplified Docker Solution
- **Final Architecture**: Single comprehensive pipeline
- **Base Image**: `continuumio/miniconda3:latest`
- **Tools**: Pre-built STAR 2.7.11b + samtools via conda
- **Scripts**: One main pipeline script (`run_full_pipeline.sh`)

### Key Technical Insights

#### Apple Silicon Compatibility
- **Issue**: Many bioinformatics tools have ARM64 compatibility problems
- **Manifestation**: Tools run but fail silently or produce 0 results
- **Docker Solution**: Provides consistent x86_64 Linux environment

#### Conda Environment Interference
- **Discovery**: Even homebrew-installed STAR affected by conda base environment
- **Root Cause**: Conda modifies PATH and library loading
- **Solution**: Docker provides isolated environment

#### Memory Management
- **Mouse Genome**: ~2.7GB, standard indexing needs 30GB+ RAM
- **Optimization**: Reduced to 15GB through sparse arrays and smaller indices
- **Trade-off**: Minimal performance impact, significant memory savings

## Project Outcome

Successfully created a robust, reproducible Docker-based pipeline that:
- ✅ Solves Apple Silicon compatibility issues
- ✅ Optimizes memory usage for common hardware
- ✅ Provides comprehensive documentation
- ✅ Maintains clean repository structure
- ✅ Enables efficient transposable element analysis

---

**Generated**: October 31, 2025
**Repository**: https://github.com/jcbclffrd/RNAseqSTAR
