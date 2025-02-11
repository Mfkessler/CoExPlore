# Persistent package, mount app directory as volume

# Use slim Python 3.13.2 as the base image
FROM python:3.13.2-slim

# Set the working directory
WORKDIR /CoExPlore

# Copy environment files
COPY requirements.txt .

# Copy package temporarily for installation
COPY package /tmp/package/

# Install system dependencies including git, graphviz and build-essential (for gcc)
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    graphviz \
    build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Pip packages
RUN pip install --no-cache-dir -r requirements.txt

# Install package
RUN pip install --no-deps /tmp/package

# Install custom packages: PyWGCNA
RUN pip install git+https://github.com/Mfkessler/PyWGCNA.git

# Install Snakemake and specific Pulp version
RUN pip install snakemake pulp==2.7.0

# Clean up temporary files and pip cache
RUN rm -rf /tmp/package && rm -rf ~/.cache/pip /tmp/* /var/tmp/*