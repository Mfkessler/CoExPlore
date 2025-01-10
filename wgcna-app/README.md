# WGCNA-App

The **WGCNA-App** is a Dockerized web application for Weighted Gene Co-Expression Network Analysis (WGCNA). It supports RNA-seq data exploration, co-expression module visualization, and integration of optional metadata (e.g., GO terms, orthogroups).

---

## Features

- Interactive exploration of WGCNA results
- Support for `.h5ad` and `.h5` TOM files
- Scalable via Celery for intensive tasks
- Ready-to-deploy with Docker Compose

---

## Quick Start

1. **Clone the repository**:
   ```bash
   git clone https://git.computational.bio.uni-giessen.de/mkessler/wgcna-app.git
   ```

2. **Build and start the app**:
   ```bash
   docker build -t ranomics_base:latest .
   docker compose -f wgcna-app/docker/docker-compose.yml --project-name wgcna_dev up
   ```

3. **Access the app**:  
   Open [http://localhost:5000](http://localhost:5000) in your browser.

---

## Preprocessing (TODO)

A Jupyter Notebook for RNA-seq data preprocessing to `.h5ad` and `.h5` formats is planned. Stay tuned!

---

## Developer

**Micha Frederick Keßler**  
**Email**: micha.kessler@computational.bio.uni-giessen.de  
**Affiliation**: Bioinformatics and Systems Biology, Justus Liebig University, Gießen  

Repository: [WGCNA-App](https://git.computational.bio.uni-giessen.de/mkessler/wgcna-app)

---

## License

Licensed under the MIT License.