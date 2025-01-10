"""
perform_WGCNA.py
-----------------
This script creates a WGCNA object for each dataset and saves the results in the output folder.

The steps include:
1. Load the count matrix and metadata
2. Create a WGCNA object
3. Perform preprocessing
4. Find modules
5. Analyze WGCNA object
6. Save WGCNA results

The script can be run in two ways:
1. Automatically: The script will automatically set the parameters based on the dataset name.
   Naming convention: {species}_{tissue}_isoform.TMM.EXPR.matrix
2. Manually: The user can set the parameters manually. 
   Example: python create_WGCNA_obj.py ../count_matrices2/AC.isoform.TMM.EXPR.matrix AC ../output/AC/
"""

import sys
import PyWGCNA
import ranomics.utils as rutils
import matplotlib.pyplot as plt
from memory_profiler import profile

@profile
def main():
    # Set matplotlib parameters
    plt.rcParams['savefig.bbox'] = 'tight'

    # Define count matrices path
    count_matrices = "/vol/share/ranomics_app/count_matrices2"
    use_isoforms = True

    # Define datasets
    names = ['PS', 'CS', 'EC', 'EG', 'HC', 'AC', 'SP', 'TT']

    # Set save options
    save_tom = True
    save_adjacency_matrix = False
    save_WGCNA = False
    figure_type = 'png'

    # Reduce the count matrix size?
    reduce_matrix = False

    # Print all datasets
    print("Starting WGCNA for datasets:", names)
    for name in names:

        """
        Prepare WGCNA object
        """

        print("Preparing WGCNA object for", name, "dataset...")
        # When using parameters
        if len(sys.argv) == 4:
            input_matrix = sys.argv[1]
            name = sys.argv[2]
            output_path = sys.argv[3] + f"{name}/"

        # Automatically set parameters
        else:
            if name == 'EC':
                input_matrix = f"{count_matrices}/{name}.gene.TMM.EXPR.matrix.updated" if not use_isoforms else f"{count_matrices}/{name}.isoform.TMM.EXPR.matrix.updated"
            else:
                input_matrix = f"{count_matrices}/{name}.gene.TMM.EXPR.matrix" if not use_isoforms else f"{count_matrices}/{name}.isoform.TMM.EXPR.matrix"

            output_path = f"/vol/share/ranomics_app/output/{name}/"

        # Print general settings
        print("Input matrix:", input_matrix)
        print("Name:", name)
        print("Output path:", output_path)

        # Print use isoforms option
        print("Use isoforms:", use_isoforms)

        # Print save options
        print("Save TOM:", save_tom)
        print("Save adjacency matrix:", save_adjacency_matrix)
        print("Save WGCNA results:", save_WGCNA)

        # Print reduce matrix option
        print("Reduce matrix:", reduce_matrix)

        # Species name
        species = rutils.find_species_by_initials(name.split("_")[0])

        # Input files: GAF and orthogroups
        if use_isoforms:
            gaf_path = f"/vol/share/ranomics_app/gaf_files/{name}.annotation.goslim_plant.gaf.gz"
            ortho_file = "/vol/share/ranomics_app/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
        else:
            gaf_path = None
            ortho_file = None

        # Output files: matrix and metadata
        matrix_file = f"{output_path}{name}_matrix.csv"
        metadata_file = f"{output_path}{name}_metadata"

        # Transform and save count table
        rutils.create_dir(output_path)
        count_df = rutils.transform_count_matrix(input_matrix)
        if reduce_matrix:
            count_df = rutils.remove_random_columns(count_df, percentage=0.9)
        count_df.to_csv(matrix_file, index=False)

        # Create and save metadata table
        metadata_df = rutils.map_tissue_types(count_df)
        metadata_df.to_csv(metadata_file, index=False)

        # Create pyWGCNA object using the matrix and metadata files
        print("Creating pyWGCNA object...")
        pyWGCNA_obj = PyWGCNA.WGCNA(name=name, 
                                    species=species, 
                                    geneExpPath=matrix_file, 
                                    outputPath=output_path,
                                    figureType=figure_type,
                                    save=True)
        # Update sample information
        pyWGCNA_obj.updateSampleInfo(path=metadata_file, sep=',')

        """
        Start WGCNA analysis
        """

        print("Starting WGCNA analysis...")

        # Perform preprocessing
        print("Preprocessing data...")
        pyWGCNA_obj.preprocess()

        # Try to find modules
        print("Finding modules...")
        pyWGCNA_obj.findModules()
        pyWGCNA_obj.datExpr.var

        # Set colors for metadata
        pyWGCNA_obj.datExpr.obs["tissue"].unique()
        pyWGCNA_obj.setMetadataColor("tissue", rutils.generate_stage_color_dict())

        print("Analyzing WGCNA object...")
        pyWGCNA_obj.analyseWGCNA()

        # Prepare and save WGCNA
        print("Preparing and saving WGCNA results...")
        rutils.prepare_and_save_wgcna(pyWGCNA_obj, output_path, gaf_path, ortho_file, save_tom=save_tom, save_adjacency_matrix=save_adjacency_matrix, save_WGCNA=save_WGCNA)
            
        print("Finished WGCNA for", name, "dataset...")

    print("Finished WGCNA for all datasets.")

if __name__ == '__main__':
    main()
