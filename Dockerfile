# Use Miniconda as the base image
FROM continuumio/miniconda3:latest

# Set the working directory
WORKDIR /CoExPlore

# Copy environment files
COPY app/requirements.txt /CoExPlore/
COPY app/environment.yml /CoExPlore/

# Copy package temporarily for installation
COPY package /tmp/package/

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    graphviz && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Conda dependencies
RUN conda env create -f environment.yml && \
    conda clean -afy

# Set the environment path
ENV PATH=/opt/conda/envs/WGCNA_APP/bin:$PATH

# Install Pip packages
RUN pip install --no-cache-dir -r requirements.txt

# Install package
RUN pip install --no-deps /tmp/package

# Install custom packages: PyWGCNA
RUN pip install git+https://github.com/Mfkessler/PyWGCNA.git

# Install Snakemake and specific Pulp version
RUN pip install snakemake pulp==2.7.0

# Clean up temporary files
RUN rm -rf /tmp/package
RUN rm -rf ~/.cache/pip /tmp/* /var/tmp/*
