# Persistent app und package

# Use Miniconda as the base image
FROM continuumio/miniconda3:latest

# Create a non-root user
RUN useradd -m appuser

# Set the working directory
WORKDIR /wgcna-app

# Copy environment files
COPY app/requirements.txt /wgcna-app/
COPY app/environment.yml /wgcna-app/

# Copy the application code
COPY app /wgcna-app/

# Copy wgcna-package temporarily for installation
COPY package /tmp/wgcna-package/

# Install Graphviz
RUN apt-get update && apt-get install -y graphviz

# Change ownership to appuser
RUN chown -R appuser:appuser /wgcna-app

# Install Conda dependencies
RUN conda env create -f environment.yml && \
    conda clean -afy

# Set the environment path
ENV PATH=/opt/conda/envs/WGCNA_APP/bin:$PATH

# Install Pip packages
RUN pip install --no-cache-dir -r requirements.txt

# Install wgcna-package
RUN pip install /tmp/wgcna-package

# Install custom packages: PyWGCNA
RUN pip install git+https://github.com/Mfkessler/PyWGCNA.git

# Clean up
RUN rm -rf /tmp/wgcna-package
RUN rm -rf ~/.cache/pip
RUN rm -rf /tmp/* /var/cache/apk/* /root/.cache/
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Switch to non-root user
USER appuser
