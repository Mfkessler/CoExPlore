# Dockerfile.app for CoExPlore services -> Hosting specific Datasets
# Persistent app and persistent package

FROM coexplore-base:latest

# Set the working directory
WORKDIR /wgcna-app

# Copy the application code
COPY app/ /wgcna-app/

# Copy the package temporarily for installation
COPY package /tmp/wgcna-package/

# Install the persistent package (without dependencies)
RUN pip install --no-deps /tmp/wgcna-package

# Clean up temporary files and pip cache
RUN rm -rf /tmp/wgcna-package ~/.cache/pip /tmp/* /var/tmp/* && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

