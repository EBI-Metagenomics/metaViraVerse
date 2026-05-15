# containers/foldseek/Dockerfile
FROM ubuntu:22.04

LABEL maintainer="metaViraVerse contributors"
LABEL description="Foldseek structural search + HHsuite for ECOD annotation"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    wget curl python3.10 python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Install Foldseek via static binary
RUN wget -q https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz && \
    tar -xzf foldseek-linux-avx2.tar.gz && \
    mv foldseek/bin/foldseek /usr/local/bin/ && \
    rm -rf foldseek foldseek-linux-avx2.tar.gz

# Install HHsuite
RUN wget -q https://github.com/soedinglab/hh-suite/releases/download/v3.3.0/hhsuite-3.3.0-SSE2-Linux.tar.gz && \
    tar -xzf hhsuite-3.3.0-SSE2-Linux.tar.gz -C /usr/local && \
    rm hhsuite-3.3.0-SSE2-Linux.tar.gz

ENV PATH="/usr/local/bin:${PATH}"

RUN pip3 install --no-cache-dir requests biopython

WORKDIR /workspace
