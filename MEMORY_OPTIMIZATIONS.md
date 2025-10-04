# Memory-Optimized STAR Parameters

## Summary of Changes

The pipeline has been updated with memory-optimized STAR parameters to work efficiently on systems with 16GB RAM.

## STAR Genome Index Building (Step 3)

**Memory-Optimized Parameters:**
- `--genomeSAsparseD 3` - Uses sparser suffix array (default is 1, higher = less RAM)
- `--genomeSAindexNbases 12` - Reduces index size (default is 14 for mouse genome)
- `--limitGenomeGenerateRAM 15000000000` - Limits RAM usage to 15GB during index generation

**Why these work:**
- `genomeSAsparseD 3` reduces the density of the suffix array, saving ~2-3GB RAM
- `genomeSAindexNbases 12` creates a smaller index suitable for the mouse genome
- `limitGenomeGenerateRAM` prevents STAR from using more than 15GB, leaving room for system

## STAR Alignment (Step 4)

**Memory-Optimized Parameters:**
- `--limitBAMsortRAM 4000000000` - Limits BAM sorting to 4GB (instead of dynamic calculation)

**Optimized for 52bp reads:**
- `--sjdbOverhang 51` - Read length - 1 (52bp - 1 = 51)
- `--outFilterMultimapNmax 100` - Allows up to 100 multimapping locations
- `--winAnchorMultimapNmax 100` - Window anchor multimap max
- `--outFilterMismatchNmax 10` - Max 10 mismatches per read
- `--outFilterScoreMinOverLread 0.3` - Min alignment score (30% of read length)
- `--outFilterMatchNmin 15` - Min 15bp matches required

## Expected Resource Usage

**Genome Index Building:**
- Time: 30-60 minutes
- Peak RAM: ~15GB
- Disk space: ~8GB for index

**Alignment:**
- Time: 30-90 minutes (depends on 135M reads)
- Peak RAM: ~10GB (genome in RAM + 4GB for sorting)
- Disk space: ~15GB for BAM file

## Verification

After building the Docker image and running the pipeline, you should see:

1. **During genome indexing:**
   ```
   Building STAR index with memory-optimized parameters...
     --genomeSAsparseD 3 (reduces RAM usage)
     --genomeSAindexNbases 12 (reduces index size)
     --limitGenomeGenerateRAM 15000000000 (15GB limit)
   ```

2. **During alignment:**
   ```
   Running STAR with memory-optimized parameters for 52bp reads...
     --limitBAMsortRAM 4000000000 (4GB for BAM sorting)
   ```

3. **Expected alignment statistics:**
   - Number of input reads: ~135,000,000
   - Uniquely mapped reads: 50-70% (typical for scRNA-seq)
   - Multi-mapped reads: 20-30%
   - Unmapped reads: <10%

## Files Updated

1. **Dockerfile** - Added comment about memory optimization
2. **scripts/run_full_pipeline.sh** - Updated STAR commands with memory-optimized parameters

## Running the Pipeline

```bash
cd /Users/jacobclifford/Desktop/scTE/mouse_analysis

# Build Docker image
./docker_build.sh

# Run pipeline
./docker_run_pipeline.sh
```

## Troubleshooting

If you still encounter memory issues:

1. **Reduce threads:** Change `THREADS=8` to `THREADS=4` in run_full_pipeline.sh
2. **Further sparse index:** Change `--genomeSAsparseD 3` to `--genomeSAsparseD 4`
3. **Smaller BAM sort:** Change `--limitBAMsortRAM 4000000000` to `--limitBAMsortRAM 2000000000`

## References

- STAR manual: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
- Memory recommendations: Section 2.2.5 of STAR manual
- Recommended for genomes >2GB with 16GB RAM
