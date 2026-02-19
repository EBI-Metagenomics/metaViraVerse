FROM mambaorg/micromamba:1.5.8

WORKDIR /opt/metacerberus

# IMPORTANT: set the target environment
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV MAMBA_DEFAULT_ENV=base
ENV PATH=/opt/conda/bin:$PATH

# Install MetaCerberus into base env
RUN micromamba install -y \
    -n base \
    -c conda-forge \
    -c bioconda \
    metacerberus \
 && micromamba clean -a -y

# Run setup (downloads databases)
RUN metacerberus.py --setup

# Install Google Chrome
RUN pip install kaleido plotly \
 && yes | plotly_get_chrome || true

CMD ["metacerberus.py", "--help"]
