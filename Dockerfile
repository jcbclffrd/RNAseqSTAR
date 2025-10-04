# Ultra-simplified Dockerfile for scTE Pipeline
# Uses conda to install STAR and samtools (avoids compilation)
# Includes memory-optimized STAR parameters for systems with 16GB RAM

FROM continuumio/miniconda3:latest

# Install STAR and samtools via conda (pre-built binaries)
RUN conda install -c bioconda -c conda-forge -y \
    star=2.7.11b \
    samtools \
    sra-tools \
    wget \
    pigz \
    && conda clean -a -y

# Install scTE for transposable element quantification
RUN pip install git+https://github.com/JiekaiLab/scTE.git

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
