FROM mambaorg/micromamba:1.5.8

SHELL ["/bin/bash", "-c"]

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

RUN pip install --no-deps bacphlip

ENTRYPOINT ["bacphlip"]