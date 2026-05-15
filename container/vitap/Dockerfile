FROM --platform=linux/amd64 condaforge/miniforge3:latest

WORKDIR /opt/vitap

# Download VITAP v1.7 from release
RUN wget https://github.com/DrKaiyangZheng/VITAP/archive/refs/tags/v.1.7.tar.gz && \
    tar -xzvf v.1.7.tar.gz && \
    mv VITAP-v.1.7 VITAP && \
    rm v.1.7.tar.gz && \
    ls -la VITAP/

# Modify the conda environment file to use diamond 0.9.14 and remove defaults channel
RUN cd VITAP && \
    sed -i 's/diamond.*/diamond=0.9.14/' VITAP_conda_environment_Linux.yaml && \
    sed -i '/- defaults/d' VITAP_conda_environment_Linux.yaml && \
    cat VITAP_conda_environment_Linux.yaml

# Remove all install_requires from setup.py (dependencies are managed by conda)
RUN cd VITAP && \
    sed -i '/install_requires/,/\]/d' setup.py && \
    echo "Modified setup.py (removed install_requires):" && \
    cat setup.py

# Create conda environment from modified yaml (without defaults channel to avoid timeout)
RUN cd VITAP && \
    conda config --remove channels defaults || true && \
    conda env create -f VITAP_conda_environment_Linux.yaml && \
    echo "conda activate VITAP" >> ~/.bashrc

# Install VITAP from source (skip dependencies - already installed via conda)
RUN cd VITAP && \
    /bin/bash -c "source activate vitap && python setup.py develop --no-deps"

# Create symlink for prodigal (conda has prodigal-gv, but VITAP expects prodigal)
RUN /bin/bash -c "source activate vitap && \
    if [ -f /opt/conda/envs/vitap/bin/prodigal-gv ]; then \
        ln -sf /opt/conda/envs/vitap/bin/prodigal-gv /opt/conda/envs/vitap/bin/prodigal; \
    fi && \
    which prodigal && prodigal -v"

# Make vitap accessible in PATH
ENV PATH="/opt/conda/envs/vitap/bin:${PATH}"

# Verify installations
RUN /bin/bash -c "source activate vitap && which VITAP && VITAP --help 2>&1 | head -5 && diamond --version"

# Create entrypoint script to auto-activate conda environment
RUN echo '#!/bin/bash' > /entrypoint.sh && \
    echo 'set -e' >> /entrypoint.sh && \
    echo 'source /opt/conda/etc/profile.d/conda.sh' >> /entrypoint.sh && \
    echo 'conda activate vitap' >> /entrypoint.sh && \
    echo 'exec "$@"' >> /entrypoint.sh && \
    chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
CMD ["/bin/bash"]
