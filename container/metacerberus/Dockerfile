FROM mambaorg/micromamba:1.5.8

USER root

# Install Chromium dependencies required by Kaleido / Chrome
RUN apt-get update && apt-get install -y \
    libnss3 \
    libatk-bridge2.0-0 \
    libcups2 \
    libxcomposite1 \
    libxdamage1 \
    libxfixes3 \
    libxrandr2 \
    libgbm1 \
    libxkbcommon0 \
    libpango-1.0-0 \
    libcairo2 \
    libasound2 \
    fonts-liberation \
    wget \
    --no-install-recommends \
 && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y \
    chromium \
    --no-install-recommends \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/metacerberus

ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV MAMBA_DEFAULT_ENV=base
ENV PATH=/opt/conda/bin:$PATH
ENV CHROME_PATH=/usr/bin/chromium

# Install MetaCerberus
RUN micromamba install -y \
    -n base \
    -c conda-forge \
    -c bioconda \
    metacerberus \
 && micromamba clean -a -y

# Setup databases
RUN metacerberus.py --setup

# Install plotly + kaleido
RUN pip install kaleido plotly

CMD ["metacerberus.py", "--help"]
