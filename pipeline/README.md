# WGCNA-Pipeline

The **WGCNA-Pipeline** is a Snakemake-based workflow for automating the preprocessing and analysis of RNA-seq data with Weighted Gene Co-Expression Network Analysis (WGCNA). 

---

## Features

- Automated pipeline for WGCNA preprocessing and analysis
- Support for optional metadata integration (e.g., GO terms, InterPro data)
- Modular Snakemake rules for flexible customization
- Docker-ready for seamless deployment

---

## Quick Start (PREVIEW!)

1. **Clone the repository**:
   ```bash
   git clone https://git.computational.bio.uni-giessen.de/mkessler/wgcna-pipeline.git
   cd wgcna-pipeline
   ```

2. **Set up the Conda environment**:
   ```bash
   conda env create -f environment.yml
   conda activate wgcna-pipeline
   ```

3. **Configure the workflow**:
   Edit the `config.yml` file to set paths and parameters for your data.

4. **Run the pipeline**:
   ```bash
   snakemake --cores all
   ```

5. **Explore the results**:  
   Check the `output/` directory for the processed data and WGCNA results.

---

## Directory Structure
- **`Snakefile`**: Defines the WGCNA preprocessing and analysis workflow.
- **`config.yml`**: Configuration file to customize the pipeline.
- **`environment.yml`**: Conda environment for dependencies.
- **`output/`**: Directory for generated results.

---

## Developer

**Micha Frederick Keßler**  
**Email**: micha.kessler@computational.bio.uni-giessen.de  
**Affiliation**: Bioinformatics and Systems Biology, Justus Liebig University, Gießen  

Repository: [WGCNA-Pipeline](https://git.computational.bio.uni-giessen.de/mkessler/wgcna-pipeline)

---

## License

Licensed under the MIT License.