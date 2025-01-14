# WGCNA-Package

The **WGCNA-Package** is a Python package designed for RNA-seq data analysis, focusing on Weighted Gene Co-Expression Network Analysis (WGCNA). It provides core functions, scripts, and tools to streamline analysis pipelines and exploratory data analysis.

---

## Features

- Core functionality for WGCNA workflows
- Jupyter notebooks for data exploration and testing
- Shell scripts for pipeline automation
- Unit tests for package validation

---

## Quick Start

1. **Clone the repository**:
   ```bash
   git clone git@git.computational.bio.uni-giessen.de:mkessler/wgcna-package.git
   cd wgcna-package
   ```

2. **Set up the Conda environment**:
   ```bash
   conda env create -f environment.yml
   conda activate wgcna-package
   bash postBuild
   ```

3. **Add the environment to Jupyter**:
   ```bash
   python -m ipykernel install --user --name=wgcna-package
   ```

4. **Explore the notebooks**:  
   Use the `notebooks/` directory to get started with data analysis and testing.

---

## Directory Structure

- **`bash/`**: Shell scripts for analysis pipelines.
- **`notebooks/`**: Jupyter notebooks for exploratory analysis.
- **`wgcna-package/`**: Core Python package.
- **`scripts/`**: Custom scripts for data processing and analysis.
- **`tests/`**: Unit tests for validating the package.

---

## Developer

**Micha Frederick Keßler**  
**Email**: micha.kessler@computational.bio.uni-giessen.de  
**Affiliation**: Bioinformatics and Systems Biology, Justus Liebig University, Gießen  

Repository: [WGCNA-Package](https://git.computational.bio.uni-giessen.de/mkessler/wgcna-package)

---

## License

Licensed under the MIT License.