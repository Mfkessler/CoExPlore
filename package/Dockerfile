# Use Miniconda as the base image
FROM continuumio/miniconda3:latest

# Set the working directory
WORKDIR /app

# Copy wgcna-package temporarily for installation
COPY . .

# Install Graphviz
RUN apt-get update && apt-get install -y graphviz

# Install Conda dependencies
RUN conda env create -f environment.yml && \
    conda clean -afy

# Set the environment path
ENV PATH=/opt/conda/envs/wgcna-package/bin:$PATH

# Install Pytest
RUN pip install pytest

# Install wgcna-package
RUN pip install .

# Install custom packages: PyWGCNA
RUN pip install git+https://github.com/Mfkessler/PyWGCNA.git

# Clean up
RUN rm -rf ~/.cache/pip
RUN rm -rf /var/cache/apk/* /root/.cache/
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the default command to run tests
CMD ["pytest", "tests/", "--disable-warnings", "--maxfail=5"]
