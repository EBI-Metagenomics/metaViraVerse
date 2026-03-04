FROM python:3.10-slim

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Full native build stack for parasail
RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    cmake \
    autoconf \
    automake \
    libtool \
    pkg-config \
    m4 \
    python3-dev \
    libc-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# Clone phamclust
RUN git clone https://github.com/chg60/phamclust.git
WORKDIR /opt/phamclust

# Upgrade build tooling
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Install deps (parasail builds here)
RUN pip install --no-cache-dir -r requirements.txt

# Install phamclust
RUN pip install --no-cache-dir .

ENTRYPOINT ["phamclust"]
