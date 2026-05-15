# containers/esmfold/Dockerfile
# ESMFold local inference container (~16 GB VRAM required for production use)
# Base: NVIDIA CUDA 11.8 + Python 3.10
FROM nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu22.04

LABEL maintainer="metaViraVerse contributors"
LABEL description="ESMFold local inference for protein structure prediction"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    python3.10 python3.10-pip python3.10-dev \
    git wget curl && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir \
    torch==2.0.1+cu118 \
    --index-url https://download.pytorch.org/whl/cu118

RUN pip3 install --no-cache-dir \
    transformers==4.36.0 \
    accelerate==0.24.0 \
    biopython==1.83 \
    requests==2.31.0

# Set HuggingFace cache to /tmp to avoid home dir issues in containers
ENV HF_HOME=/tmp/hf_cache

WORKDIR /workspace
