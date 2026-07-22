FROM mambaorg/micromamba:1.5.8

SHELL ["/bin/bash", "-c"]

# Install system tools required by Nextflow processes
USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        original-awk \
        dateutils \
        procps \
        pcregrep \
        sed \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create conda environment
RUN micromamba create -y -n bacphlip \
    -c conda-forge \
    -c bioconda \
    python=3.8 \
    numpy=1.19.5 \
    scikit-learn=0.23.1 \
    biopython=1.77 \
    pandas \
    joblib \
    hmmer \
    pip \
    && micromamba clean --all --yes

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/envs/bacphlip/bin:$PATH

# Install BACPHLIP
RUN micromamba run -n bacphlip pip install --no-deps bacphlip

# Important: remove ENTRYPOINT for Nextflow compatibility
# Nextflow overrides CMD but not ENTRYPOINT safely
CMD ["bacphlip"]
