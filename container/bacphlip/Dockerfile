FROM mambaorg/micromamba:1.5.8

SHELL ["/bin/bash", "-c"]

# Create BACPHLIP environment with explicit channels
RUN micromamba create -y -n bacphlip \
    -c conda-forge \
    -c bioconda \
    -c defaults \
    python=3.8 \
    numpy=1.19.5 \
    scikit-learn=0.23.1 \
    biopython=1.77 \
    pandas \
    joblib \
    hmmer \
    pip \
    && micromamba clean --all --yes

# Auto-activate environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/envs/bacphlip/bin:$PATH

# Install BACPHLIP
RUN pip install bacphlip

# Sanity check
RUN python - <<EOF
import numpy, sklearn, bacphlip
print("NumPy:", numpy.__version__)
print("sklearn:", sklearn.__version__)
print("BACPHLIP import OK")
EOF

ENTRYPOINT ["bacphlip"]
