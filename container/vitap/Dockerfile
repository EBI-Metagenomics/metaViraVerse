FROM condaforge/miniforge3:latest

WORKDIR /opt/vitap

# Install vitap but ignore diamond dependency upgrade
RUN conda install -c bioconda -c conda-forge vitap -y

# Download and install diamond v0.9.14 manually
RUN wget https://github.com/bbuchfink/diamond/releases/download/v0.9.14/diamond-linux64.tar.gz && \
    tar -xvf diamond-linux64.tar.gz && \
    mv diamond /usr/local/bin/diamond && \
    chmod +x /usr/local/bin/diamond && \
    rm -rf diamond-linux64.tar.gz

# Confirm version
RUN diamond --version

CMD ["/bin/bash"]
