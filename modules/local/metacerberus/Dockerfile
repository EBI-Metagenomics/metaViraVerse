FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /opt

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    bash \
    && rm -rf /var/lib/apt/lists/*

# Install mamba into base environment (required by install script)
RUN conda install -y -c conda-forge mamba

# Clone MetaCerberus
RUN git clone https://github.com/raw-lab/MetaCerberus.git

# Install MetaCerberus
WORKDIR /opt/MetaCerberus
RUN bash install_metacerberus.sh

# Make MetaCerberus conda environment default in Docker
SHELL ["conda", "run", "-n", "MetaCerberus", "/bin/bash", "-c"]
ENV PATH=/opt/conda/envs/MetaCerberus/bin:$PATH

# Default command
CMD ["metacerberus.py", "--help"]
