import os
import random
import h5py
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import networkx as nx
import igraph as ig
import re
from .ortho import add_ortho_count, update_ortho_ID
from matplotlib import cm, colors
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.godag_plot import plot_gos
from contextlib import redirect_stdout, redirect_stderr
from typing import List, Dict, Tuple, Union, Callable
from anndata import AnnData
from sklearn.impute import KNNImputer
from sklearn.preprocessing import scale
from pandas.api.types import is_numeric_dtype
from jinja2 import Environment, FileSystemLoader
from collections import defaultdict
from scipy.spatial.distance import squareform


def create_dir(path: str) -> None:
    """
    Creates a directory at the specified path if it does not already exist.

    Parameters:
    - path (str): The file path of the directory to be created.
    """

    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Directory created at: {path}")
    else:
        print(f"Directory already exists at: {path}")


def find_species_by_initials(species: str, file_path: str = "../info/species") -> str:
    """
    Finds and returns the full name of a species based on the initial letters of its name provided in a single string.

    Parameters:
    - species (str): A string containing the initial letters of the first and second word of the species name.
    - file_path (str): The path to the file containing the species names.

    Returns:
    - str: The full name of the species that matches the given initials. If no match is found, returns None.
    """

    if len(species) != 2:
        return species

    first_letter, second_letter = species[0].upper(), species[1].upper()

    with open(file_path, 'r') as file:
        for line in file:
            words = line.strip().split()
            if len(words) >= 2 and words[0][0].upper() == first_letter and words[1][0].upper() == second_letter:
                return line.strip()

    raise ValueError(
        f"No species found with the provided initials. Please check the species file: {os.path.abspath(file_path)}")


def remove_column_if_exists(df: pd.DataFrame, column_name: str) -> None:
    """
    Removes a specified column from a pandas DataFrame if it exists.

    Parameters:
    - df (pd.DataFrame): DF from which the column will be removed.
    - column_name (str): Name of the column to be removed.
    """

    if column_name in df.columns:
        df.drop(columns=[column_name], inplace=True)
        print(f"Column '{column_name}' was removed.")
    else:
        print(f"Column '{column_name}' not found.")


def transform_count_matrix(file_path: str) -> pd.DataFrame:
    """
    Transforms a count matrix from a format where transcripts are rows and samples are columns
    to a format where each row represents a sample and columns represent transcripts.

    The input file should be a tab-separated values (TSV) file with transcripts as rows and samples as columns.
    The first row should contain sample names, and the first column should contain identifiers.

    The transformed data frame is saved as a CSV file at the specified output path without the index.

    Parameters:
    - file_path (str): The path to the input TSV file containing the count matrix.

    Returns:
    - pd.DataFrame: Table with samples as rows and transcripts as columns.
    """

    df = pd.read_csv(file_path, sep='\t')
    transformed_df = df.T
    transformed_df.columns = transformed_df.iloc[0]
    transformed_df = transformed_df.drop(transformed_df.index[0])
    transformed_df.reset_index(inplace=True)
    transformed_df.rename(columns={'index': 'sample'}, inplace=True)

    return transformed_df


def remove_random_columns(df: pd.DataFrame, percentage: float = 0.9, seed: int = 0) -> pd.DataFrame:
    """
    Removes a random percentage of columns from a pandas DataFrame, excluding the first column.

    Parameters:
    - df (pd.DataFrame): DataFrame from which columns will be removed.
    - percentage (float): The fraction of columns to remove, between 0 and 1.
    - seed (int): Seed for the random number generator for reproducibility.

    Returns:
    - pd.DataFrame: Table with the specified percentage of columns removed.
    """

    np.random.seed(seed)
    num_cols_to_remove = int(np.ceil((len(df.columns) - 1) * percentage))
    cols_to_remove = np.random.choice(
        df.columns[1:], size=num_cols_to_remove, replace=False)

    return df.drop(columns=cols_to_remove)


def get_tissue_map(metadata: str = "RanOmics") -> Dict[str, str]:
    """
    Returns a dictionary mapping tissue codes to their corresponding tissue types.

    Returns:
    - Dict[str, str]: A dictionary where keys are tissue codes and values are tissue types.
    """

    if metadata == "RanOmics":
        tissue_map = {
            "A": "Bud stage 1",
            "B": "Bud stage 2",
            "C": "Bud stage 3",
            "D": "Bud stage 4",
            "E": "Sepals at anthesis",
            "F": "Petals at anthesis",
            "G": "Stamens at anthesis",
            "H": "Gynoecia at anthesis",
            "I": "Shoot apex",
            "J": "Petal early stage",
            "K": "Petal mid stage",
            "L": "Petal late stage",
            "M": "Young fruits",
            "N": "Mid-stage fruit",
            "O": "Seeds 5 dap",
            "P": "Root",
            "Q": "Young leaf",
            "R": "Mature leaf",
            "S": "Seedling",
            "T": "Mature petal nectary part",
            "U": "Mature petal no nectary part",
            "V": "Non-spurred Petal"
        }
    elif metadata == "Ceratopteris":
        tissue_map = {
            # TODO: Add tissue mapping for Ceratopteris
        }

    return tissue_map


def map_tissue_types(df: pd.DataFrame, sample_col: str = "sample", tissue_col: str = "tissue") -> pd.DataFrame:
    """
    This function reads a specified column from the input DataFrame and creates a new DataFrame with only
    the specified column and a new 'tissue' column. It maps substrings found between two "_" characters in the
    specified column to their corresponding tissue types according to a predefined mapping table.

    Parameters:
    - df (pd.DataFrame): The input DataFrame.
    - sample_col (str): The name of the column to read from the input DataFrame.
    - tissue_col (str): The name of the new column to be added to the output DataFrame.

    Returns:
    - pd.DataFrame: A new table with only the specified and the mapped 'tissue' column.
    """

    tissue_map = get_tissue_map()

    def map_sample_to_tissue(sample):
        for key in tissue_map.keys():
            if f"_{key}_" in sample:
                return tissue_map[key]
        return "Unknown"

    mapped_tissues = df[sample_col].apply(map_sample_to_tissue)
    new_df = pd.DataFrame(
        {sample_col: df[sample_col], tissue_col: mapped_tissues})

    return new_df


def get_general_info(adata: AnnData, obs_col: str = "tissue") -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    This function calculates the number of samples, genes, and tissues in an AnnData object.

    Parameters:
    - adata (AnnData): The AnnData object to analyze.
    - obs_col (str): The name of the column in adata.obs containing tissue information.

    Returns:
    - pd.DataFrame: A DataFrame containing the number of samples, genes, tissues, and other info.
    """

    species = adata.uns['species']
    num_samples = adata.n_obs
    num_transcripts = adata.n_vars
    num_modules = adata.var['moduleColors'].nunique()
    num_tissues = adata.obs[obs_col].nunique()
    total_counts = adata.X.sum()

    # numbers of samples per tissue
    samples_per_tissue = adata.obs[obs_col].value_counts().to_dict()

    data = {
        "Species": [species],
        "Samples": [num_samples],
        "Transcripts": [num_transcripts],
        "Total Counts": [round(total_counts)],
        "Modules": [num_modules],
        "Tissues": [num_tissues],
    }

    df_info = pd.DataFrame(data)
    df_samples_per_tissue = pd.DataFrame(
        list(samples_per_tissue.items()), columns=["Tissue", "Samples"])

    return df_info, df_samples_per_tissue


def get_info_for_list(adatas: List[AnnData], obs_col: str = "tissue") -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    This function applies get_general_info to each AnnData object in a list.

    Parameters:
    - adatas (List[AnnData]): The list of AnnData objects to analyze.
    - obs_col (str): The name of the column in adata.obs containing tissue information.

    Returns:
    - pd.DataFrame: A DataFrame containing the combined info for each AnnData object.
    """

    all_info = []
    all_samples_per_tissue = []

    for adata in adatas:
        temp_obs_col = obs_col if obs_col in adata.obs else "Combined_Trait"

        if temp_obs_col not in adata.obs:
            raise ValueError(f"{temp_obs_col} not found in adata.obs")

        df_info, df_samples_per_tissue = get_general_info(adata, temp_obs_col)
        all_info.append(df_info)
        all_samples_per_tissue.append(df_samples_per_tissue)

    combined_info = pd.concat(all_info).reset_index(drop=True)
    combined_samples_per_tissue = pd.concat(
        all_samples_per_tissue).reset_index(drop=True)

    return combined_info, combined_samples_per_tissue


def compute_sample_statistics_df(adata: AnnData) -> pd.DataFrame:
    """
    Calculates statistical measures for each sample in the count matrix of an AnnData object
    and returns them in a pandas DataFrame.

    Parameters:
    - adata (AnnData): An AnnData object containing the count matrix in adata.X.

    Returns:
    - pd.DataFrame: Table containing the calculated statistical measures,
      with samples as rows and statistical measures as columns.
    """

    X = adata.X

    # Check if X is a sparse matrix and convert to a dense matrix if needed
    if sp.issparse(X):
        X_dense = X.toarray()
    else:
        X_dense = X

    # Calculate statistical measures across rows (samples)
    stats = {
        'Mean': np.mean(X_dense, axis=1),
        'Median': np.median(X_dense, axis=1),
        'Standard Deviation': np.std(X_dense, axis=1),
        'Max': np.max(X_dense, axis=1),
        'Min': np.min(X_dense, axis=1),
        'Sum': np.sum(X_dense, axis=1)
    }

    # Convert the stats dictionary to a DataFrame
    stats_df = pd.DataFrame(stats, index=adata.obs_names)

    return stats_df


def generate_stage_color_dict(custom_stages: List[str] = None) -> Dict[str, str]:
    """
    Generates a dictionary for setting metadata colors, based on predefined stages.
    Includes up to 23 distinct stages plus an 'Unknown' category, using a combination
    of matplotlib's tab20 and Set3 colormaps for distinct colors.

    Parameters:
    - custom_stages (List[str]): A list of custom stages to use instead of the predefined stages.

    Returns:
    - Dict[str, str]: Dictionary where keys are the stage names and values are distinct colors.
    """

    # Predefined stages for Ranomics
    if not custom_stages:
        stages = [
            "Bud stage 1", "Bud stage 2", "Bud stage 3", "Bud stage 4",
            "Sepals at anthesis", "Petals at anthesis", "Stamens at anthesis",
            "Gynoecia at anthesis", "Shoot apex", "Petal early stage",
            "Petal mid stage", "Petal late stage", "Young fruits",
            "Mid-stage fruit", "Seeds 5 dap", "Root", "Young leaf",
            "Mature leaf", "Seedling", "Mature petal nectary part",
            "Mature petal no nectary part", "Non-spurred Petal", "Unknown"
        ]
    # Custom stages
    else:
        stages = custom_stages

    if len(stages) > 30:
        raise ValueError(
            "This function supports up to 30 distinct categories.")

    # Generate colors from matplotlib colormaps
    # First 20 colors from tab20
    colors_tab20 = cm.tab20(np.linspace(0, 1, 20))
    colors_set3 = cm.Set3(np.linspace(0, 1, 12))    # First 12 colors from Set3

    # Combine the two colormaps and take only as many colors as needed
    all_colors = np.vstack((colors_tab20, colors_set3))[:len(stages)]

    # Convert the RGB colors to hex
    hex_colors = [colors.to_hex(color) for color in all_colors]

    # Create the color dictionary
    color_dict = dict(zip(stages, hex_colors))

    return color_dict


def process_anndata(
    adata: AnnData,
    gaf_gz_path: str = None,
    ortho_dir: str = None,
    ortho_column_name: str = None,
    decimal_places: int = 2,
    go_terms: bool = True,
    ortho_id: bool = True,
    qc_metrics: bool = True
) -> None:
    """
    Enhances an AnnData object with additional genomic annotations and recalculates QC metrics.
    If decimal_places is specified, the float columns in the AnnData object are rounded to the specified number of decimal places.

    This function performs the following operations:
    1. Adds gene ontology (GO) terms to the AnnData object.
    2. Adds orthologous IDs to the AnnData object.
    3. Calculates quality control metrics for the object.
    4. Drops the 'n_cells_by_counts' column.
    5. Rounds float columns to the specified number of decimal places. (Optional)

    Parameters:
    - adata (AnnData): The AnnData object to be processed.
    - gaf_gz_path (str): Path to the Gene Annotation File (GAF) in gzipped format.
    - ortho_dir (str): Path to the dir containing orthogroups.
    - ortho_column_name (str): Column name under which orthologous IDs are stored. E.g. species abbreviation.
    - decimal_places (int): Number of decimal places to round float columns to.
    - go_terms (bool): If True, adds GO terms to the AnnData object.
    - ortho_id (bool): If True, adds orthologous IDs to the AnnData object.
    - qc_metrics (bool): If True, recalculates quality control metrics for the AnnData object.
    """

    if go_terms and gaf_gz_path:
        add_go_terms_to_anndata(adata, gaf_gz_path=gaf_gz_path)
    if ortho_id and ortho_column_name:
        add_ortho_ids_to_anndata(
            adata, column_name=ortho_column_name, base_path=ortho_dir)
        update_ortho_ID(adata)
    if qc_metrics:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        # Drop n_cells_by_counts column
        adata.var.drop(columns="n_cells_by_counts", axis=1, inplace=True)
    if decimal_places:
        float_cols = adata.var.select_dtypes(include=['float64']).columns
        adata.var[float_cols] = adata.var[float_cols].round(decimal_places)


def add_go_terms_to_anndata(adata: AnnData, gaf_gz_path: str, hog: bool = False) -> None:
    """
    Reads a GAF.gz (Gene Annotation File) and adds the GO (Gene Ontology) terms 
    associated with each gene to the .var attribute of an AnnData object. The 
    function assumes that the AnnData object's .var attribute contains gene 
    identifiers as the index.

    The GAF file is expected to be compressed in gzip format. The function 
    extracts the GO terms and the corresponding gene symbols from the GAF file, 
    then merges this information into the .var table of the specified AnnData object
    as a comma-separated string.

    Parameters:
    - adata (AnnData): An AnnData object containing gene expression data.
    - gaf_gz_path (str): Path to the GAF.gz file containing gene annotations and GO terms.
    - hog (bool): If True, the function assumes that the identifiers in the .var are HOG IDs.
    """

    if hog:
        mapping_df = pd.read_csv(gaf_gz_path, sep='\t', header=None, names=[
                                 'hog_id', 'go_terms'])
        hog_to_go_dict = pd.Series(
            mapping_df.go_terms.values, index=mapping_df.hog_id).to_dict()
        adata.var['go_terms'] = adata.var.index.map(hog_to_go_dict).fillna('')
    else:
        # Read GAF file and extract necessary columns
        gaf_df = pd.read_csv(gaf_gz_path, sep='\t', header=None, comment='!', compression='gzip',
                             usecols=[1, 4], names=['Transcript', 'GO_ID'])

        # Group by transcript and aggregate GO terms into a comma-separated string
        gene_go_df = gaf_df.groupby('Transcript')['GO_ID'].agg(
            lambda x: ','.join(x)).reset_index()
        # Rename the 'GO_ID' column to 'go_terms'
        gene_go_df.columns = ['Transcript', 'go_terms']

        # Merge the GO terms into the '.var' table of AnnData
        adata.var = adata.var.merge(
            gene_go_df, left_index=True, right_on='Transcript', how='left').set_index('Transcript')

        # Count the number of GO terms for each transcript
        adata.var['n_go_terms'] = adata.var['go_terms'].apply(
            lambda x: 0 if pd.isna(x) else x.count(',') + 1)


def aggregate_go_terms(series: pd.Series) -> str:
    """
    Aggregate GO terms, remove duplicates, and return them as a comma-separated string.

    Parameters:
    - pd.Series: A pandas Series containing GO terms as a comma-separated string.

    Returns:
    - str: A comma-separated string containing unique GO terms.
    """

    terms = series.dropna().astype(str).str.split(',')
    all_terms = [term.strip()
                 for sublist in terms for term in sublist if term.strip()]
    unique_terms = np.unique(all_terms)

    return ','.join(unique_terms) if unique_terms.size > 0 else ""


def filter_adata_by_obs(adata: AnnData, obs_column: str = "tissue", expr_threshold: float = None, display_stats: bool = True) -> AnnData:
    """
    Filters an AnnData object by applying a variable expression threshold for each unique value
    in a specified .obs column, setting the counts of transcripts that don't meet the threshold to zero,
    and finally removes transcripts that are zero across all samples. Additionally, reports the number of transcripts
    excluded per tissue and overall, both as absolute numbers and percentages.

    Parameters:
    - adata (AnnData): The AnnData object containing .var, .obs, and .X (expression matrix).
    - obs_column (str): The name of the .obs column to analyze by (e.g., 'tissue').
    - expr_threshold (float): A fixed threshold for minimum expression. If None, the threshold
                                        is calculated for each unique value in the obs_column.
    - display_stats (bool): If True, display the stats of excluded transcripts.

    Returns:
    - AnnData: The modified AnnData object where non-expressed transcripts are set to zero and entirely non-expressed transcripts are removed.
    """

    if obs_column not in adata.obs:
        obs_column = "Combined_Trait"
        if obs_column not in adata.obs:
            raise ValueError(f"{obs_column} not found in adata.obs")

    unique_values = adata.obs[obs_column].unique()
    exclusion_counts = np.zeros((adata.n_vars, len(unique_values)), dtype=bool)

    # Apply thresholds and update expression
    for i, value in enumerate(unique_values):
        tissue_mask = adata.obs[obs_column] == value
        filtered_adata = adata[tissue_mask, :]

        current_expr_threshold = calculate_expr_threshold(
            filtered_adata) if expr_threshold is None else expr_threshold
        is_expressed = filtered_adata.X.sum(axis=0) > current_expr_threshold

        # Update exclusion counts and set non-expressed counts to zero
        exclusion_counts[:, i] = ~is_expressed
        adata.X[tissue_mask, :][:, ~is_expressed] = 0

    # Check for transcripts that are zero across all samples
    all_zero_mask = np.all(exclusion_counts, axis=1)
    adata = adata[:, ~all_zero_mask]

    # Create statistics for excluded transcripts
    exclusion_stats = {
        value: {
            'Excluded transcripts': np.sum(exclusion_counts[:, i]),
            'Percentage excluded': np.mean(exclusion_counts[:, i]) * 100
        } for i, value in enumerate(unique_values)
    }
    exclusion_stats['Overall'] = {
        'Excluded transcripts': np.sum(all_zero_mask),
        'Percentage excluded': np.mean(all_zero_mask) * 100
    }

    # Transpose for better readability
    exclusion_df = pd.DataFrame(exclusion_stats).T

    if display_stats:
        print(exclusion_df)

    return adata


def calculate_expr_threshold(adata: AnnData, std_multiplier: float = 1) -> float:
    """
    Calculate the expression threshold for 'safely expressed' transcripts based on the median and standard deviation of expression values.

    Parameters:
    - adata (AnnData): The AnnData object containing expression data.
    - std_multiplier (float): The multiplier for the standard deviation to adjust the threshold.

    Returns:
    - float: The calculated expression threshold.
    """

    # Ensure that the data is in a dense format for computation
    if not isinstance(adata.X, np.ndarray):
        X_dense = adata.X.toarray()
    else:
        X_dense = adata.X

    # Calculate the median and standard deviation of the expression values across all cells
    median_expression = np.median(X_dense, axis=0)
    std_expression = np.std(X_dense, axis=0)

    # Calculate the threshold for 'safely expressed'
    threshold = np.median(median_expression) + \
        std_multiplier * np.std(std_expression)

    return threshold


def add_goea_to_anndata(adata: AnnData, obo_path: str = "../go-basic.obo", results_dir: str = "goea_results", go_column: str = "go_terms",
                        uns_key: str = "goea_results", column: str = 'moduleColors', top_percentage: int = 10, go_background_set: set = None) -> None:
    """
    Performs GOEA using GO terms from an AnnData object, integrating detailed results and
    visualizations directly into the AnnData object. Considers only the top percentage of 
    highest expressed genes per module based on average expression across samples.

    Parameters:
    - adata (AnnData): The AnnData object containing transcripts and their associated GO terms.
    - obo_path (str): Path to the .obo file containing the GO ontology.
    - results_dir (str): Directory to save GOEA visualizations and results.
    - go_column (str): Column in adata.var containing comma separated GO terms for each transcript.
    - uns_key (str): Key in AnnData.uns to store GOEA results. If None, results are not saved in AnnData.uns.
    - column (str): Column in adata.var containing module labels or other groupings for GOEA.
    - top_percentage (int): Percentage of top expressed genes to consider in each module.
    - go_background_set (set): Set of all genes to consider as background for GOEA. If None, all genes in adata.var are used.
    """

    if top_percentage < 1:
        top_percentage = 1
        print("Top percentage set to 1 to ensure at least one transcript is selected.")

    # Prepare a dummy file to capture all outputs
    with open(os.devnull, 'w') as fnull:
        with redirect_stdout(fnull), redirect_stderr(fnull):
            go_dag = GODag(obo_path)
            os.makedirs(results_dir, exist_ok=True)

            # Prepare the GO background set
            if go_background_set is None:
                go_background_set = set(adata.var.index)

            # Create a dictionary mapping transcript IDs to GO IDs
            assoc = {}
            for transcript_id, go_ids_str in adata.var[go_column].dropna().items():
                go_ids = go_ids_str.split(",")
                assoc[transcript_id] = set(go_ids)

            # Initialize GOEA object
            goea_obj = GOEnrichmentStudyNS(
                pop=go_background_set,
                ns2assoc={"generic_namespace": assoc},
                godag=go_dag,
                methods=['fdr_bh'])

            adata.uns[uns_key] = {}

            # Iterate over modules
            for module_label in adata.var[column].unique():
                module_transcripts = adata.var.loc[adata.var[column]
                                                   == module_label]

                # Select top percentage of transcripts based on mean expression
                top_cutoff = int(len(module_transcripts)
                                 * (top_percentage / 100))
                top_transcripts = module_transcripts.nlargest(
                    top_cutoff, 'mean_counts')

                # Get all transcript IDs for the current module after filtering
                module_transcript_ids = top_transcripts.index.tolist()

                if not module_transcript_ids:
                    continue

                # Perform GOEA
                goea_results = goea_obj.run_study(module_transcript_ids)

                # Filter significant results and extract additional information
                goea_results_sig = [(r.GO, r.name, r.p_fdr_bh, r.ratio_in_study[0]/r.ratio_in_study[1],
                                     r.ratio_in_study[0], go_dag[r.GO].depth) for r in goea_results if r.p_fdr_bh < 0.05]

                # Plot and save results if there are significant findings
                if goea_results_sig:
                    plot_path = os.path.join(
                        results_dir, f'goea_plot_{module_label}.png')
                    significant_go_terms = [r[0] for r in goea_results_sig]
                    plot_gos(plot_path, significant_go_terms, go_dag,
                             title=f"{adata.uns['name']}: GO Terms for {module_label}")

                    adata.uns[uns_key][module_label] = {
                        'significant_terms': goea_results_sig,
                        'plot_path': plot_path
                    }


def get_goea_df(adata: AnnData, key: str = 'goea_results', column_name: str = "Tissue") -> pd.DataFrame:
    """
    Extracts GOEA results from an AnnData object and returns a DataFrame with detailed results.

    Parameters:
    - adata (AnnData): The AnnData object containing GOEA results in adata.uns[key].
    - key (str): The key in adata.uns where the GOEA results are stored.
    - column_name (str): The name of the column to use for the module labels in the resulting DataFrame.

    Returns:
    - pd.DataFrame: A DataFrame with detailed GOEA results.
    """

    goea_results = adata.uns[key]
    goea_results_list = []

    for module_label, data in goea_results.items():
        for term in data['significant_terms']:
            go_id, term_name, p_value, fold_enrichment, num_genes, depth = term
            goea_results_list.append({
                column_name: module_label,
                "GO_ID": go_id,
                "Term": term_name,
                "P-Value": p_value,
                "Fold Enrichment": fold_enrichment,
                "Num Genes": num_genes,
                "Depth": depth,
                # Add plot path if needed
                "Plot Path": data.get('plot_path', 'N/A')
            })

    return pd.DataFrame(goea_results_list)


def add_ortho_id_to_anndata(adata: AnnData, mapping_file: str, column_name: str) -> None:
    """
    Adds a new column 'ortho_ID' to the adata.var dataframe by mapping identifiers
    from a specified column of a given mapping file directly to the adata.var of an AnnData object.

    Parameters:
    - adata (AnnData): The AnnData object whose .var dataframe will be modified.
    - mapping_file (str): Path to the mapping file containing orthogroup IDs and identifiers.
    - column_name (str): Name of the column in the mapping file that contains the identifiers to be mapped.
    """

    if not os.path.exists(mapping_file):
        print(
            f"Warning: Mapping file '{mapping_file}' not found. Filling 'ortho_ID' with empty strings.")
        adata.var['ortho_ID'] = ''
        return

    ortho_df = pd.read_csv(mapping_file, sep="\t", usecols=[
                           "HOG", column_name], index_col="HOG")

    adata.var['ortho_ID'] = None

    id_to_ortho = {}
    for ortho_id, identifiers in ortho_df[column_name].dropna().items():
        for identifier in identifiers.split(", "):
            id_to_ortho[identifier] = ortho_id

    adata.var['ortho_ID'] = adata.var.index.map(id_to_ortho).fillna(None)
    add_ortho_count(adata)


def add_ortho_ids_to_anndata(adata: AnnData, column_name: str, base_path: str = "../Phylogenetic_Hierarchical_Orthogroups") -> None:
    """
    Adds new columns 'N0_ortho', 'N1_ortho', ..., 'N6_ortho' to the adata.var dataframe by mapping identifiers
    from specified columns of the given mapping files. Additionally, sets the 'ortho_ID' column to the values of 'N0_ortho'.

    Parameters:
    - adata (AnnData): The AnnData object whose .var dataframe will be modified.
    - column_name (str): Name of the column in the mapping files that contains the identifiers to be mapped.
    - base_path (str): Base path to the directory containing orthogroup mapping files (N0.tsv, N1.tsv, etc.).
    """

    for level in range(7):
        node_key = f"N{level}_ortho"
        mapping_file = f"{base_path}/N{level}.tsv"
        ortho_df = pd.read_csv(mapping_file, sep="\t", usecols=[
                               "HOG", column_name], index_col="HOG")

        id_to_ortho = {}
        for ortho_id, identifiers in ortho_df[column_name].dropna().items():
            for identifier in identifiers.split(", "):
                if identifier not in id_to_ortho:
                    id_to_ortho[identifier] = []
                id_to_ortho[identifier].append(ortho_id)

        adata.var[node_key] = adata.var.index.map(
            lambda x: ','.join(id_to_ortho.get(x, []))).fillna("")

    adata.var['ortho_ID'] = adata.var['N0_ortho']
    add_ortho_count(adata)


def analyze_adata_by_obs_and_module(adata: AnnData, obs_column: str = None) -> Dict[str, pd.DataFrame]:
    """
    Analyze AnnData by module for each unique value in a specified .obs column or for the entire adata if obs_column is None.

    Parameters:
    - adata (AnnData): AnnData containing .var, .obs, and .X (expression matrix) attributes.
    - obs_column (str): The name of the .obs column to analyze by (e.g., 'tissue').

    Returns:
    - A dictionary of DataFrames, each key is a unique value from the obs_column or 'all_data' if obs_column is None,
      and each value is a DataFrame with 'moduleColors' as the index and two columns:
      'transcripts_per_module' and 'unique_GO_terms_per_module'.
    """

    results = {}
    adata.var['GO_ID'] = adata.var['GO_ID'].apply(
        lambda x: eval(x) if isinstance(x, str) else x)

    def calculate_module_stats(adata: AnnData) -> pd.DataFrame:
        return adata.var.groupby('moduleColors').agg(
            transcripts_per_module=('moduleColors', 'size'),
            unique_GO_terms_per_module=('GO_ID', lambda x: len(
                set().union(*[set(item) for item in x.dropna()])))
        ).reset_index()

    if obs_column and obs_column in adata.obs:
        unique_values = adata.obs[obs_column].unique()
        for value in unique_values:
            # Filter AnnData object for each unique value in .obs column and calculate module stats
            filtered_adata = adata[adata.obs[obs_column] == value].copy()
            is_expressed = np.asarray(filtered_adata.X.sum(
                axis=0) > calculate_expr_threshold(filtered_adata)).flatten()
            filtered_adata = filtered_adata[:, is_expressed].copy()
            results[value] = calculate_module_stats(filtered_adata)

    else:
        # Filter and Calculate statistics for the entire adata
        is_expressed = np.asarray(adata.X.sum(
            axis=0) > calculate_expr_threshold(adata)).flatten()
        adata = adata[:, is_expressed].copy()
        results['all_data'] = calculate_module_stats(adata)

    return results


def filter_single_sample_groups(adata: AnnData, obs_column: str = "tissue") -> AnnData:
    """
    Filters out groups from an anndata object that have only one sample in the specified
    observation (obs) column.

    Parameters:
    - adata (AnnData): The anndata object to filter.
    - obs_column (str): The name of the column in adata.obs based on which to filter groups.

    Returns:
    - AnnData: The filtered anndata object with groups having more than one sample.
    """

    if obs_column not in adata.obs:
        obs_column = "Combined_Trait"
        if obs_column not in adata.obs:
            raise ValueError(f"{obs_column} not found in adata.obs")

    groups_with_single_sample = adata.obs[obs_column].value_counts() == 1
    groups_with_single_sample = groups_with_single_sample[groups_with_single_sample].index.tolist(
    )

    indices_to_remove = adata.obs.index[adata.obs[obs_column].isin(
        groups_with_single_sample)]
    adata_filtered = adata[~adata.obs.index.isin(indices_to_remove)].copy()

    # Print the number of removed samples and their indices
    print(
        f"Removed {len(indices_to_remove)} samples: {indices_to_remove.tolist()}")

    return adata_filtered


def update_adata_with_dge_info(adata: AnnData, rank_df: pd.DataFrame = None, uns_group_key: str = "DGE", uns_score_key: str = "DGE_score") -> None:
    """
    Updates an anndata object with differential gene expression (DGE) group information and scores,
    based on a ranking DataFrame. The function sorts the DataFrame, removes duplicates based on names,
    and maps the unique groups and scores back to the anndata object's .var table.

    Parameters:
    - adata (AnnData): The anndata object to be updated.
    - rank_df (pd.DataFrame): A DataFrame containing at least "names", "group", and "scores"
      columns, where "names" correspond to variable names in adata.var.
    """

    if rank_df is None:
        rank_df = sc.get.rank_genes_groups_df(adata, group=None)

    # Sort and deduplicate the rank DataFrame
    rank_df_sorted = rank_df.sort_values(
        by=["names", "scores"], ascending=[True, False])
    rank_df_dedup = rank_df_sorted.drop_duplicates(
        subset=["names"], keep="first")

    # Map groups and scores to the var names
    names_to_group = pd.Series(
        rank_df_dedup.group.values, index=rank_df_dedup.names).to_dict()
    names_to_scores = pd.Series(
        rank_df_dedup.scores.values, index=rank_df_dedup.names).to_dict()

    # Update the anndata var DataFrame with DGE group and score
    adata.var[uns_group_key] = adata.var.index.map(
        names_to_group).fillna('No group')
    adata.var[uns_score_key] = adata.var.index.map(names_to_scores).fillna(0)


def get_top_expressed_genes_by_tissue(adatas: List[AnnData], n: int = 1, column: str = None, trait: str = "tissue") -> List[Dict[str, List[str]]]:
    """
    Generates a list of dictionaries, each containing the top N expressed genes for each tissue from an AnnData object.
    The keys of each dictionary are in the format <Count>_<Tissue>_<Species>.

    Parameters
    ----------
    - adatas (List[AnnData]): List of AnnData objects containing expression data.
    - n (int): Number of top expressed genes to display per tissue.
    - column (str): The name of the column in adata.var to return as gene identifiers, or None to return the index.
    - trait (str): The name of the column in adata.obs containing the tissue information.

    Returns
    -------
    - List[Dict[str, List[str]]]: A list of dictionaries where each dictionary contains tissues as keys
      in the format <Tissue>_<Species> and lists of top N expressed genes as values.
    """

    gene_dicts = []

    for adata in adatas:
        if trait not in adata.obs:
            trait = "Combined_Trait"
            if trait not in adata.obs:
                raise ValueError(f"{trait} not found in adata.obs")

        species_name = adata.uns['species']

        tissue_genes = {}
        tissues = adata.obs[trait].unique()

        for tissue in tissues:
            tissue_data = adata[adata.obs[trait] == tissue].copy()
            # Sort genes based on mean expression and select the top N
            top_gene_indices = np.argsort(tissue_data.X.mean(axis=0))[-n:]
            if column:
                # Retrieve the desired gene identifiers using the column parameter
                top_genes = adata.var[column].iloc[top_gene_indices].tolist()
            else:
                # Retrieve the index if column is None
                top_genes = adata.var.iloc[top_gene_indices].index.tolist()

            key = f"{tissue}_{species_name}"
            tissue_genes[key] = top_genes

        gene_dicts.append(tissue_genes)

    return gene_dicts


def create_mapping_dict(adata: AnnData, mapping_column: str = "moduleColors", value_column: str = None,
                        transcripts: Dict[str, List[str]] = None, min_count: int = None) -> dict:
    """
    Creates a dictionary mapping entries from the specified mapping column in adata.var to a list of associated transcripts
    or values from another specified column. Only the specified transcripts for each species are considered.

    If mapping_column is set to "species", the dictionary will aggregate all data in adata.var without considering
    a specific column, and the key will be <count>_<species>.

    Parameters
    ----------
    - adata (AnnData): An AnnData object containing transcript data.
    - mapping_column (str): The column in adata.var to use for creating the mapping. If set to "species", the entire dataset is considered.
    - value_column (str): The column from which values are extracted to form the list associated with each key in the mapping. 
      If None, the transcript identifiers are used.
    - transcripts (Dict[str, List[str]]): A dictionary with species as keys and lists of transcripts to consider. 
      If None, all transcripts are considered.
    - min_count (int): Minimum number of values associated with a key to be included in the mapping.

    Returns
    -------
    - dict: A dictionary where keys are unique entries from the mapping_column and values are lists of unique, non-empty associated 
      transcripts or values from the value_column.
    """

    mapping_dict = {}

    # Filter based on species-specific transcripts if provided
    species = adata.uns['species']
    if transcripts is not None and species in transcripts:
        adata_var_filtered = adata.var.loc[adata.var.index.isin(
            transcripts[species])]
    else:
        adata_var_filtered = adata.var

    if adata_var_filtered.empty:
        print("Warning: No matching transcripts found after filtering. Returning an empty dictionary.")
        return mapping_dict

    if mapping_column == "species":
        species_name = adata.uns['species']
        values = []

        for transcript, row in adata_var_filtered.iterrows():
            if value_column:
                if pd.notna(row[value_column]):
                    value = row[value_column].strip()
                    if value:  # Filter out empty strings
                        values.append(value)
            else:
                values.append(transcript)

        # Remove duplicates and empty values
        values = list(set(filter(None, values)))
        count = len(values)
        key = f"{count}_{species_name}"
        mapping_dict[key] = values

    else:
        for transcript, row in adata_var_filtered.iterrows():
            if pd.notna(row[mapping_column]):
                elements = row[mapping_column].split(',')
                for element in elements:
                    # Include species in the element key
                    element = f"{element}_{species}"
                    if element not in mapping_dict:
                        mapping_dict[element] = []

                    if value_column:
                        if pd.notna(row[value_column]):
                            value = row[value_column].strip()
                            if value:  # Filter out empty strings
                                mapping_dict[element].append(value)
                    else:
                        mapping_dict[element].append(transcript)

        # Remove duplicates and empty values, and adjust key names with counts
        for key in list(mapping_dict.keys()):
            mapping_dict[key] = list(set(filter(None, mapping_dict[key])))

            count = len(mapping_dict[key])
            new_key = f"{count}_{key}"
            mapping_dict[new_key] = mapping_dict.pop(key)

    # Filter dictionary by min_count if specified
    if min_count:
        mapping_dict = {key: value for key,
                        value in mapping_dict.items() if len(value) >= min_count}

    return mapping_dict


def calculate_jaccard_index(set1: set, set2: set) -> float:
    """
    Calculates the Jaccard index between two sets.

    Parameters:
    - set1 (set): The first set.
    - set2 (set): The second set.

    Returns:
    - float: The Jaccard index between the two sets.
    """

    if not set1 or not set2:
        return 0.0

    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))

    return intersection / union if union != 0 else 0.0


def get_jaccard_df(dicts: List[Dict[str, List[str]]], suffixes: List[str] = None) -> pd.DataFrame:
    """
    Creates a DataFrame containing Jaccard indices from a list of dictionaries.

    Parameters
    ----------
    - dicts (List[Dict[str, List[str]]]): List of dictionaries.
    - suffixes (List[str], optional): List of suffixes corresponding to each dictionary. If not provided,
      the suffixes will be extracted from the keys in the dictionaries.

    Returns
    -------
    - pd.DataFrame: A DataFrame containing Jaccard indices between each pair.
    """

    # If suffixes are not provided, extract them from the keys in the dictionaries
    if suffixes is None:
        suffixes = []
        for dic in dicts:
            # Extract the suffix from the first key
            suffix = next(iter(dic)).split('_')[-1]
            suffixes.append(suffix)

    if len(dicts) != len(suffixes):
        raise ValueError(
            "The number of dictionaries must match the number of suffixes.")

    # Create a combined dictionary with module colors as keys
    module_colors = {}
    for dic in dicts:
        for key in dic.keys():
            module_colors[key] = dic[key]

    # Flatten the list of module colors for DataFrame indexing
    flat_module_colors = list(module_colors.keys())

    # Initialize an empty DataFrame for Jaccard indices
    jaccard_df = pd.DataFrame(
        0, index=flat_module_colors, columns=flat_module_colors, dtype=float)

    # Calculate Jaccard index for each pair of module colors from each dictionary combination
    for color1 in flat_module_colors:
        set1 = set(module_colors[color1])
        for color2 in flat_module_colors:
            set2 = set(module_colors[color2])
            jaccard_index = calculate_jaccard_index(set1, set2)
            jaccard_df.at[color1, color2] = jaccard_index

    return jaccard_df


def extract_top_similarities(jaccard_df: pd.DataFrame, top_n: int = 5, threshold: float = 0.4, inter_species_only: bool = False) -> pd.DataFrame:
    """
    Extracts the top N highest similarities from a Jaccard similarity DataFrame, with optional filtering to exclude intra-species comparisons.

    Parameters:
    - jaccard_df (pd.DataFrame): The input DataFrame containing the Jaccard similarity indices.
    - top_n (int): Number of top similarities to return.
    - threshold (float): Minimum similarity value to consider.
    - inter_species_only (bool): If True, only inter-species comparisons are considered.

    Returns:
    - pd.DataFrame: A DataFrame of the top N similarities, excluding intra-species comparisons if specified.
    """

    # Create a mask to identify values above the threshold and off the diagonal
    mask = np.triu(np.ones_like(jaccard_df, dtype=bool), k=1)
    mask &= jaccard_df > threshold

    # Apply mask and stack the DataFrame to get a Series with MultiIndex
    filtered_series = jaccard_df.where(mask).stack()

    if inter_species_only:
        # Extract species from row and column names
        row_species = filtered_series.index.get_level_values(
            0).map(lambda x: x.split('_')[-1])
        col_species = filtered_series.index.get_level_values(
            1).map(lambda x: x.split('_')[-1])

        # Filter out intra-species comparisons
        filtered_series = filtered_series[row_species != col_species]

    top_values = filtered_series.sort_values(ascending=False).head(top_n)

    result_df = pd.DataFrame(top_values, columns=['Jaccard Value'])
    result_df.index.names = ['Row Name', 'Column Name']

    return result_df.reset_index()


def prepare_and_save_wgcna(pyWGCNA_obj: object, output_path: str, gaf_path: str = None, ortho_file: str = None,
                           save_tom: bool = False, save_adjacency_matrix: bool = False, save_WGCNA: bool = False,
                           n: bool = 2, qc_metrics: bool = True) -> None:
    """
    Prepare and save the WGCNA data analysis results.

    This function processes and saves various outputs from a WGCNA object including gene expressions, 
    eigengenes, hub genes, and matrices like TOM and adjacency. It also handles file writing for 
    both the main data and orthologous data after processing.

    Parameters:
    - pyWGCNA_obj (object): The pyWGCNA object containing WGCNA analysis results and configurations.
    - output_path (str): The base path where all the result files will be saved.
    - gaf_path (str): Path to the GAF file containing gene annotations and GO terms.
    - ortho_file (str): Path to the file containing orthologous gene IDs.
    - save_tom (bool): If True, saves the topological overlap matrix.
    - save_adjacency_matrix (bool): If True, saves the adjacency matrix.
    - save_WGCNA (bool): If True, saves the whole WGCNA object.
    - create_ortho (bool): If True, creates an orthologous AnnData object and saves it.
    - n (int): Number of decimal places to round float columns to.
    - qc_metrics (bool): If True, calculates QC metrics for the AnnData object.
    """

    if save_WGCNA:
        # Save whole WGCNA object
        pyWGCNA_obj.saveWGCNA()

    # Prepare adata from WGCNA object
    adata = pyWGCNA_obj.datExpr
    name = pyWGCNA_obj.name
    hub_genes = pyWGCNA_obj.top_n_hub_genes(n=100)
    hub_dict = {color: data for color,
                data in hub_genes.groupby('moduleColors')}
    for color, df in hub_dict.items():
        df.index = df.index.get_level_values(1)

    adata.raw = pyWGCNA_obj.geneExpr
    adata.obsm["eigengenes"] = pyWGCNA_obj.datME
    adata.varm["eigengenes"] = pyWGCNA_obj.signedKME

    # Store additional metadata
    adata.uns["name"] = name
    adata.uns["output_path"] = os.path.abspath(pyWGCNA_obj.outputPath)
    adata.uns["power"] = pyWGCNA_obj.power
    adata.uns["species"] = pyWGCNA_obj.species
    adata.uns["network_type"] = pyWGCNA_obj.networkType
    adata.uns["TOM_type"] = pyWGCNA_obj.TOMType
    adata.uns["min_module_size"] = pyWGCNA_obj.minModuleSize
    adata.uns["sft"] = pyWGCNA_obj.sft.astype(float)
    adata.uns["gene_tree"] = pyWGCNA_obj.geneTree
    adata.uns["metadata_corr"] = {
        "module_corr": pyWGCNA_obj.moduleTraitCor, "module_p": pyWGCNA_obj.moduleTraitPvalue}
    adata.uns["hub_genes"] = hub_dict
    adata.uns["num_transcripts"] = pyWGCNA_obj.geneExpr.shape[1]

    # Process adata
    if gaf_path:
        add_go_terms_to_anndata(adata, gaf_gz_path=gaf_path)
    if ortho_file:
        add_ortho_id_to_anndata(adata, ortho_file, column_name=name)
    if qc_metrics:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        # Drop n_cells_by_counts column
        adata.var.drop(columns="n_cells_by_counts", axis=1, inplace=True)

    # Round float columns to n decimal places
    float_cols = adata.var.select_dtypes(include=['float64']).columns
    adata.var[float_cols] = adata.var[float_cols].round(n)

    # Save the annotated adata and ortho data
    adata.write_h5ad(f"{output_path}{name}.h5ad")

    if save_tom:
        # Save TOM outside of adata
        print("Saving topological overlap...")
        save_matrix_to_hdf5(
            pyWGCNA_obj.TOM, f"{output_path}tom_matrix.h5", save_sym_only=False, include_diagonal=True)

    if save_adjacency_matrix:
        # Save adjacency matrix outside of adata
        print("Saving adjacency matrix...")
        save_matrix_to_hdf5(pyWGCNA_obj.adjacency,
                            f"{output_path}adjacency_matrix.h5")


def get_adatas_info_df(adatas: List[AnnData], adatas_ortho: List[AnnData], trait: str = "tissue") -> pd.DataFrame:
    """
    Extracts and combines information from multiple AnnData objects into a DataFrame.

    Parameters:
    - adatas (List[AnnData]): List of AnnData objects.
    - adatas_ortho (List[AnnData]): List of orthologous AnnData objects.
    - trait (str): The name of the column in adata.obs to use for grouping.

    Returns:
    - pd.DataFrame: A DataFrame containing information from both the original and orthologous AnnData objects.
    """

    combined_data = []

    for adata, adata_ortho in zip(adatas, adatas_ortho):
        adatas_df = pd.DataFrame({
            'Name': [adata.uns['species']],
            'Samples': [adata.n_obs],
            'Transcripts': [adata.n_vars],
            'Total Counts': [adata.obs['total_counts'].sum()],
            'Modules': [adata.var['moduleLabels'].nunique()],
            'Tissues': [adata.obs[trait].nunique()],
            'Go Terms': [adata.var['go_terms'].apply(lambda x: len(x.split(',')) if pd.notna(x) else 0).sum()],
        })

        adatas_df['Go Terms'] = adatas_df['Go Terms'].astype(int)

        adatas_ortho_df = pd.DataFrame({
            'Name': [adata_ortho.uns['species'] + " Ortho"],
            'Samples': [adata_ortho.n_obs],
            'Transcripts': [adata_ortho.n_vars],
            'Total Counts': [adata_ortho.obs['total_counts'].sum()],
            'Modules': [adata_ortho.var['moduleLabels'].nunique()],
            'Tissues': [adata_ortho.obs[trait].nunique()],
            'Go Terms': [adata_ortho.var['go_terms'].apply(lambda x: len(x.split(',')) if pd.notna(x) else 0).sum()],
        })

        adatas_ortho_df['Go Terms'] = adatas_ortho_df['Go Terms'].astype(int)

        combined_df = pd.concat([adatas_df, adatas_ortho_df])
        combined_data.append(combined_df)

    final_df = pd.concat(combined_data)

    return final_df


def reduce_tom_matrix(TOM: pd.DataFrame, reduction_percentage: float) -> pd.DataFrame:
    """
    Reduces the TOM matrix by randomly removing a given percentage of genes.

    Parameters:
    - TOM (pd.DataFrame): The TOM matrix.
    - reduction_percentage (float): The percentage of genes to remove.

    Returns:
    - pd.DataFrame: The reduced TOM matrix.
    """

    num_genes = TOM.shape[0]
    num_to_remove = int(num_genes * reduction_percentage / 100)
    genes_to_remove = random.sample(list(TOM.index), num_to_remove)

    reduced_TOM = TOM.drop(genes_to_remove, axis=0).drop(
        genes_to_remove, axis=1)

    return reduced_TOM


def save_matrix_to_hdf5(table: pd.DataFrame, filename: str, save_sym_only: bool = False, include_diagonal: bool = True) -> None:
    """
    Saves the table (pd.DataFrame) into an HDF5 file with optimized row-wise chunking.
    
    Parameters:
    - table (pd.DataFrame): The table to be saved.
    - filename (str): The name of the HDF5 file.
    - save_sym_only (bool): If True, saves only the upper triangle of a symmetric matrix.
    - include_diagonal (bool): When saving only the upper triangle, specifies whether to include
      the diagonal elements.
    
    The function saves the index and column labels, and either the full data matrix or only
    the upper triangle, depending on the value of 'save_sym_only'. In the full-matrix mode,
    the data is saved with row-wise chunking (i.e. each row is stored in a separate chunk)
    to optimize later row-based access.
    """

    with h5py.File(filename, 'w') as f:
        # Save index and column labels
        f.create_dataset('index', data=np.array(table.index, dtype='S'))
        f.create_dataset('columns', data=np.array(table.columns, dtype='S'))

        if save_sym_only:
            # Ensure the matrix is square
            n_rows, n_cols = table.shape
            assert n_rows == n_cols, "The matrix must be square to save only half."

            # Extract the upper triangle indices
            if include_diagonal:
                triu_indices = np.triu_indices(n_rows, k=0)
            else:
                triu_indices = np.triu_indices(n_rows, k=1)

            # Flatten the upper triangle data
            upper_triangle = table.values[triu_indices]

            # Save the upper triangle data
            f.create_dataset('data_upper', data=upper_triangle)
            f.attrs['n'] = n_rows
            f.attrs['save_sym_only'] = True
            f.attrs['include_diagonal'] = include_diagonal
        else:
            # Save the full data matrix with optimized row-wise chunking.
            data = table.values
            chunk_size = (1, data.shape[1])  # Jeder Chunk entspricht einer Zeile.
            f.create_dataset('data', data=data, chunks=chunk_size)
            f.attrs['save_sym_only'] = False
            f.attrs['include_diagonal'] = True


def load_subset_from_hdf5(filename: str, rows: list, cols: list = None, threshold: float = None,
                          index_map=None, columns_map=None, use_symmetry: bool = False,
                          include_neighbours: bool = False, max_neighbors: int = 10) -> Tuple[pd.DataFrame, set]:
    """
    Loads a subset of a large HDF5 file into a pandas DataFrame.
    Optionally, if include_neighbours is True, for each transcript in rows the function
    adds neighbors with connection value >= threshold (up to max_neighbors).
    
    Parameters:
    - filename (str): HDF5 file name.
    - rows (list): List of transcript labels (for both rows and columns).
    - cols (list): List of transcript labels for columns. If None, uses rows.
    - threshold (float): Only edges with value >= threshold are considered.
    - index_map (dict): Mapping from transcript labels to row indices.
    - columns_map (dict): Mapping from transcript labels to column indices.
    - use_symmetry (bool): If True, uses symmetry (not supported with include_neighbours).
    - include_neighbours (bool): If True, add neighbors with connection >= threshold.
    - max_neighbors (int): Maximum number of neighbors to add per transcript.
    
    Returns:
    - pd.DataFrame: The filtered TOM subset.
    - set: Set of neighbor transcripts (if include_neighbours is True), else empty set.
    """
    original_rows = set(rows)
    rows = list(set(rows))
    if cols:
        cols = list(set(cols))
    else:
        cols = rows

    neighbor_set = set()
    if include_neighbours:
        if threshold is None:
            raise ValueError("A threshold must be provided if include_neighbours is True.")
        if use_symmetry:
            raise NotImplementedError("include_neighbours is not supported with use_symmetry=True.")
        if index_map is None or columns_map is None:
            with h5py.File(filename, 'r') as f:
                file_index = [s.decode('utf-8') for s in f['index'][:]]
                file_columns = [s.decode('utf-8') for s in f['columns'][:]]
            index_map = {label: idx for idx, label in enumerate(file_index)}
            columns_map = {label: idx for idx, label in enumerate(file_columns)}
        new_transcripts = set(rows)
        inv_columns_map = {v: k for k, v in columns_map.items()}
        with h5py.File(filename, 'r') as f:
            data = f['data']
            for transcript in rows:
                if transcript not in index_map:
                    continue
                r_idx = index_map[transcript]
                row_data = data[r_idx, :]
                neighbor_idxs = np.where(row_data >= threshold)[0]
                neighbor_idxs = neighbor_idxs[neighbor_idxs != r_idx]
                if len(neighbor_idxs) > max_neighbors:
                    sorted_idx = neighbor_idxs[np.argsort(row_data[neighbor_idxs])[::-1]]
                    neighbor_idxs = sorted_idx[:max_neighbors]
                for col_idx in neighbor_idxs:
                    neighbor = inv_columns_map.get(col_idx)
                    if neighbor:
                        new_transcripts.add(neighbor)
        neighbor_set = new_transcripts - original_rows
        rows = list(new_transcripts)
        cols = list(new_transcripts)

    if index_map is None or columns_map is None:
        with h5py.File(filename, 'r') as f:
            file_index = [s.decode('utf-8') for s in f['index'][:]]
            file_columns = [s.decode('utf-8') for s in f['columns'][:]]
        index_map = {label: idx for idx, label in enumerate(file_index)}
        columns_map = {label: idx for idx, label in enumerate(file_columns)}

    try:
        row_indices = [index_map[row] for row in rows]
        col_indices = [columns_map[col] for col in cols]
    except KeyError as e:
        raise Exception("Label not found in index or columns: " + str(e))

    sorted_row_order = np.argsort(row_indices)
    sorted_col_order = np.argsort(col_indices)
    row_indices_sorted = np.array(row_indices)[sorted_row_order]
    col_indices_sorted = np.array(col_indices)[sorted_col_order]

    try:
        with h5py.File(filename, 'r') as f:
            if use_symmetry:
                data_upper = f['data_upper'][:]
                n = f.attrs['n']
                full_matrix = squareform(data_upper)
                if full_matrix.shape[0] < n:
                    full_matrix = np.pad(full_matrix, ((0, 1), (0, 1)),
                                         mode='constant', constant_values=0)
                block_data = full_matrix[np.ix_(row_indices_sorted, col_indices_sorted)]
            else:
                data = f['data']
                data_rows = data[row_indices_sorted, :]
                block_data = data_rows[:, col_indices_sorted]
    except Exception as e:
        raise Exception("Error loading block from HDF5 file: " + str(e))

    original_row_order = np.argsort(sorted_row_order)
    original_col_order = np.argsort(sorted_col_order)
    data_subset = block_data[original_row_order, :][:, original_col_order]

    df_subset = pd.DataFrame(data_subset, index=rows, columns=cols)

    if threshold is not None:
        np.fill_diagonal(df_subset.values, np.nan)
        mask = df_subset > threshold
        df_subset = df_subset.where(mask, other=np.nan)
        df_subset = df_subset.dropna(how='all', axis=0).dropna(how='all', axis=1)

    return df_subset, neighbor_set


def convert_hdf5_to_upper_triangle(old_filename, new_filename):
    """
    Converts a full square matrix stored in an HDF5 file into a flattened upper triangle matrix
    and saves it into a new HDF5 file.

    Parameters:
    - old_filename (str): The filename of the original HDF5 file containing the full matrix.
    - new_filename (str): The filename for the new HDF5 file to store the upper triangle matrix.
    """

    with h5py.File(old_filename, 'r') as f_old:
        data = f_old['data'][:]
        index_labels = [s.decode('utf-8') for s in f_old['index'][:]]
        column_labels = [s.decode('utf-8') for s in f_old['columns'][:]]

    n = data.shape[0]
    assert n == data.shape[1], "The matrix must be square."

    # Extract the upper triangle matrix without the diagonal
    triu_indices = np.triu_indices(n, k=1)
    upper_triangle = data[triu_indices]

    with h5py.File(new_filename, 'w') as f_new:
        f_new.create_dataset('data_upper', data=upper_triangle)
        f_new.attrs['n'] = n
        f_new.create_dataset('index', data=np.array(index_labels, dtype='S'))
        f_new.create_dataset(
            'columns', data=np.array(column_labels, dtype='S'))


def get_tom_maps(file_path: str) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Extracts index and column maps from an HDF5 file containing a TOM matrix.

    Parameters:
    - file_path (str): The path to the HDF5 file containing the TOM matrix.

    Returns:
    - Tuple[Dict[str, List[str]], Dict[str, List[str]]]: A tuple containing the row and column maps.
    """

    with h5py.File(file_path, 'r') as f:
        index = [s.decode('utf-8') for s in f['index'][:]]
        columns = [s.decode('utf-8') for s in f['columns'][:]]

    index_map = {label: idx for idx, label in enumerate(index)}
    columns_map = {label: idx for idx, label in enumerate(columns)}

    return index_map, columns_map


def add_MADS(adatas: Union[AnnData, List[AnnData]], path: str = "../info/MADS") -> None:
    """
    Add a 'MADS' column to the var attribute of AnnData objects.

    Parameters:
    - adatas (Union[AnnData, List[AnnData]]): Single AnnData object or list of AnnData objects.
    - path (str): Path to the MADS TSV file containing the plant names and transcripts.

    Returns:
    - None
    """

    mads_df = pd.read_csv(path, sep='\t', header=None,
                          names=['Plant', 'Transcript'])

    if isinstance(adatas, AnnData):
        adatas = [adatas]

    # Process each AnnData object
    for adata in adatas:
        if 'MADS' not in adata.var.columns:
            adata.var['MADS'] = False

        plant_name = adata.uns['name']
        plant_mads = mads_df[mads_df['Plant'] == plant_name]['Transcript']
        adata.var['MADS'] = adata.var.index.isin(plant_mads).astype(bool)


def load_comparison(adatas: List[AnnData], directory: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Loads the matrices and returns only the subsets that pertain to the species in the provided adatas.

    Parameters:
    - adatas (List[AnnData]): List of AnnData objects.
    - directory (str): Directory from which the files will be loaded.

    Returns:
    - Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: Filtered correlation matrix, p-value matrix, and combined DataFrame.
    """

    species_list = [adata.uns['species'] for adata in adatas]

    # Load matrices
    corr = pd.read_csv(f"{directory}/correlation_matrix.csv", index_col=0)
    pvals = pd.read_csv(f"{directory}/p_value_matrix.csv", index_col=0)
    combined_df = pd.read_csv(f"{directory}/combined_df.csv", index_col=0)

    # Filter columns based on species names
    filtered_columns_corr_pvals = [col for col in corr.columns if any(
        species in col for species in species_list)]
    filtered_columns_combined = [col for col in combined_df.columns if any(
        species in col for species in species_list)]

    # Filter rows and columns for correlation and p-value matrices
    corr_filtered = corr.loc[filtered_columns_corr_pvals,
                             filtered_columns_corr_pvals]
    pvals_filtered = pvals.loc[filtered_columns_corr_pvals,
                               filtered_columns_corr_pvals]

    # Filter columns for combined DataFrame
    combined_df_filtered = combined_df[filtered_columns_combined]

    return corr_filtered, pvals_filtered, combined_df_filtered


def identify_network_clusters(G: nx.Graph, cluster_name: str,
                              node_threshold: int = None, node_threshold_percent: float = 0.02,
                              print_info: bool = False) -> Dict[str, str]:
    """
    Assigns nodes to clusters based on the graph object and returns a mapping of nodes to cluster assignments.

    Parameters:
    - G (nx.Graph): The graph object.
    - cluster_name (str): Prefix for cluster assignment.
    - node_threshold (int, optional): Minimum number of nodes required in a cluster to be assigned.
    - node_threshold_percent (float): Percentage of nodes required in a cluster to be assigned.
    - print_info (bool): If True, prints cluster information.

    Returns:
    - Dict[str, str]: Mapping of transcript IDs to cluster assignments.
    """

    if not isinstance(cluster_name, str):
        raise TypeError("cluster_name must be a string")

    # Find connected components
    clusters = list(nx.connected_components(G))

    if not node_threshold:
        if not node_threshold_percent:
            raise ValueError(
                "Either node_threshold or node_threshold_percent must be provided")

        node_threshold = int(len(G.nodes) * node_threshold_percent)

    # Assign clusters to transcripts if they meet the node_threshold
    filtered_clusters = [
        cluster for cluster in clusters if len(cluster) >= node_threshold]

    cluster_map = {}
    for new_cluster_id, cluster in enumerate(filtered_clusters, start=1):
        for transcript in cluster:
            cluster_map[transcript] = f"{cluster_name} {new_cluster_id}"

    # Fill in transcripts that are not in any cluster
    for transcript in G.nodes:
        if transcript not in cluster_map:
            cluster_map[transcript] = "No clusters"

    if print_info:
        cluster_counts = pd.Series(cluster_map).value_counts().sort_index()
        for cluster_id, count in cluster_counts.items():
            cluster_transcripts = [
                t for t, c in cluster_map.items() if c == cluster_id]
            print(
                f"{cluster_id}: {count} transcripts - {', '.join(cluster_transcripts)}")

    return cluster_map


def add_at_id(adata: AnnData, filename: str = "/vol/ranomics/mads-box/N0.tsv") -> None:
    """
    Creates a dictionary mapping specified plant transcripts to AT gene names with transcript versions based on orthogroups from a given file,
    and adds a new column 'AT_ID' to adata.var with the AT transcript names corresponding to each specified plant transcript.

    Parameters:
    - adata (AnnData): The AnnData object containing the data.
    - filename (str): The path to the mapping file.
    """

    name = adata.uns['name']

    def parse_mapping_line(line: str) -> dict:
        parts = line.split('\t')
        orthogroup = parts[1]
        at_transcripts = []
        plant_transcripts = []

        for part in parts:
            if 'AT_' in part:
                at_transcripts.extend(part.split())
            if f'{name}_' in part:
                plant_transcripts.extend(part.split())

        return {orthogroup: {'AT': at_transcripts, name: plant_transcripts}}

    def extract_gene(at_transcript: str) -> str:
        match = re.match(r"AT_([A-Z0-9]+)", at_transcript)
        if match:
            return match.group(1)
        return at_transcript

    def remove_prefix(transcript: str) -> str:
        return re.sub(r'^[A-Z]{2}_', '', transcript)

    mapping_dict = {}

    with open(filename, 'r') as file:
        for line in file:
            if line.strip():
                parsed_line = parse_mapping_line(line)
                for orthogroup, transcripts in parsed_line.items():
                    for plant_transcript in transcripts[name]:
                        plant_transcript_no_prefix = remove_prefix(
                            plant_transcript)
                        if plant_transcript_no_prefix not in mapping_dict:
                            mapping_dict[plant_transcript_no_prefix] = []
                        gene_names = [extract_gene(remove_prefix(
                            at_transcript)) for at_transcript in transcripts['AT']]
                        mapping_dict[plant_transcript_no_prefix].extend(
                            gene_names)

    def map_transcripts(plant_transcript):
        if plant_transcript in mapping_dict:
            at_transcripts = mapping_dict[plant_transcript]
            return ','.join(at_transcripts)
        return None

    adata.var['AT_ID'] = adata.var.index.map(map_transcripts)
    adata.var['AT_ID'] = adata.var['AT_ID'].str.replace(
        r'\.Araport[^\s,]*', '', regex=True)

# Based on moduleEigengenes function from PyWGCNA "https://github.com/mortazavilab/PyWGCNA/blob/main/PyWGCNA/wgcna.py"


def calculate_eigengenes(adata: AnnData, cluster_map: Dict[str, str], column: str = "clusters", exclude: List[str] = None,
                         impute: bool = True, nPC: int = 1, align: str = "along average", scaleVar: bool = True,
                         store_in_obsm: bool = False) -> dict:
    """
    Calculates module eigengenes (1st principal component) of modules in a given single dataset.

    Parameters:
    - adata (AnnData): Annotated data matrix.
    - cluster_map (Dict[str, str]): Mapping of genes to module assignments.
    - column (str): Column name in adata.var containing module assignments (used for reference only, no changes).
    - exclude (list): List of modules to exclude from the eigengene calculation.
    - impute (bool): Whether to impute missing data.
    - nPC (int): Number of principal components to calculate.
    - align (str): Align eigengenes with average expression.
    - scaleVar (bool): Scale expression data before SVD.
    - store_in_obsm (bool): Whether to store eigengenes in adata.obsm.

    Returns:
    - dict: Dictionary containing "eigengenes", "averageExpr", "varExplained", and "moduleInfo".
    """

    expr = pd.DataFrame(adata.X, index=adata.obs_names,
                        columns=adata.var_names)
    module_colors = adata.var['moduleColors']

    # Remove alle Keys of the cluster_map that are not in the adata.var.index
    cluster_map = {key: value for key,
                   value in cluster_map.items() if key in expr.columns}

    clusters = pd.Series(expr.columns.map(cluster_map).fillna(
        "No clusters"), index=expr.columns, name=column)
    modules = pd.Categorical(clusters).categories

    if exclude:
        modules = [mod for mod in modules if mod not in exclude]

    PrinComps = pd.DataFrame(index=expr.index)
    averExpr = pd.DataFrame(index=expr.index)
    varExpl = pd.DataFrame(index=range(nPC), columns=modules)
    moduleInfo = {}

    included_genes = []

    for mod in modules:
        module_genes = expr.columns[clusters == mod]
        expr_mod = expr[module_genes]

        # if mod is not "No clusters" include the gene
        if mod != "No clusters":
            included_genes.extend(module_genes)

        mod = f"{mod}"

        # Imputation and scaling
        if impute and expr_mod.size > 0:
            imputer = KNNImputer(n_neighbors=min(10, expr_mod.shape[0] - 1))
            expr_mod = imputer.fit_transform(expr_mod)

        if scaleVar and expr_mod.size > 0:
            expr_mod = pd.DataFrame(
                scale(expr_mod, axis=0), index=expr.index, columns=module_genes)

        # Calculate SVD and align eigengene
        if expr_mod.size > 0:
            u, s, vh = np.linalg.svd(expr_mod, full_matrices=False)
            eigengene = u[:, 0]
            if align == "along average":
                avg_expr = expr_mod.mean(axis=1)
                if np.corrcoef(eigengene, avg_expr)[0, 1] < 0:
                    eigengene = -eigengene
            PrinComps[str(mod)] = eigengene
            averExpr["AE" + str(mod)] = expr_mod.mean(axis=1)
            varExpl.loc[:, mod] = (s[:nPC]**2) / np.sum(s**2)

        # Calculate module info
        total_transcripts = len(module_genes)
        color_counts = module_colors[module_genes].value_counts()
        most_common_module = color_counts.idxmax() if not color_counts.empty else None
        most_common_module_count = color_counts.max() if not color_counts.empty else 0
        unique_modules = color_counts[color_counts > 0].count()

        moduleInfo[mod] = {
            "species": adata.uns['species'],
            "total_transcripts": total_transcripts,
            "most_common_module": most_common_module,
            "most_common_module_percentage": round((most_common_module_count / total_transcripts) * 100, 2) if total_transcripts > 0 else 0,
            "unique_modules": unique_modules
        }

    # Calculate for "All Clusters"
    all_dataModule = expr[included_genes]
    total_transcripts_all = len(included_genes)

    if total_transcripts_all > 0:
        color_counts_all = module_colors[included_genes].value_counts()
        most_common_module_all = color_counts_all.idxmax(
        ) if not color_counts_all.empty else None
        most_common_module_count_all = color_counts_all.max(
        ) if not color_counts_all.empty else 0
        unique_modules_all = color_counts_all[color_counts_all > 0].count()

        moduleInfo["All Clusters"] = {
            "species": adata.uns['species'],
            "total_transcripts": total_transcripts_all,
            "most_common_module": most_common_module_all,
            "most_common_module_percentage": round((most_common_module_count_all / total_transcripts_all) * 100, 2) if total_transcripts_all > 0 else 0,
            "unique_modules": unique_modules_all
        }

        # Impute and scale if applicable
        if impute and not all_dataModule.empty:
            imputer = KNNImputer(n_neighbors=min(
                10, all_dataModule.shape[0] - 1))
            all_dataModule = imputer.fit_transform(all_dataModule)

    if scaleVar and all_dataModule.size > 0:
        all_dataModule = pd.DataFrame(
            scale(all_dataModule, axis=0), index=expr.index, columns=included_genes)

    # Calculate SVD and align eigengene for all clusters
    if all_dataModule.size > 0:
        u, s, vh = np.linalg.svd(all_dataModule, full_matrices=False)
        eigengene = u[:, 0]
        if align == "along average":
            avg_expr = all_dataModule.mean(axis=1)
            if np.corrcoef(eigengene, avg_expr)[0, 1] < 0:
                eigengene = -eigengene
        PrinComps[f"{adata.uns['name']}: All Clusters"] = eigengene
        averExpr[f"{adata.uns['name']}: AE_All"] = all_dataModule.mean(axis=1)
        varExpl[f"{adata.uns['name']}: All Clusters"] = (
            s[:nPC]**2) / np.sum(s**2)

    # Store results in adata.obsm if specified
    if store_in_obsm:
        adata.obsm[column] = PrinComps

    # Convert moduleInfo to DataFrame
    moduleInfo_df = pd.DataFrame(moduleInfo).T
    moduleInfo_df.index.name = "Cluster"

    # Remove the "No clusters" key if it exists
    PrinComps = PrinComps.drop(f"No clusters", axis=1, errors='ignore')
    moduleInfo_df = moduleInfo_df.drop(f"No clusters", axis=0, errors='ignore')

    return {
        "eigengenes": PrinComps,
        "averageExpr": averExpr,
        "varExplained": varExpl,
        "moduleInfo": moduleInfo_df
    }

# Based on getDatTraits function from PyWGCNA "https://github.com/mortazavilab/PyWGCNA/blob/main/PyWGCNA/wgcna.py"


def get_data_trait_df(adata: AnnData, meta_data: str) -> pd.DataFrame:
    """
    Converts metadata from an AnnData object to a DataFrame with binary traits.

    Parameters:
    - adata (AnnData): The AnnData object containing metadata.
    - meta_data (str): The name of the metadata column to convert.

    Returns:
    - pd.DataFrame: A DataFrame with binary traits for each unique value in the metadata column.
    """

    data = adata.obs.copy()[meta_data]
    datTraits = pd.DataFrame(index=data.index)
    for i in range(data.shape[1]):
        if is_numeric_dtype(data.iloc[:, i].dtypes):
            datTraits[data.columns[i]] = data.iloc[:, i]
            continue
        data.iloc[:, i] = data.iloc[:, i].astype(str)
        if len(np.unique(data.iloc[:, i])) == 2:
            datTraits[data.columns[i]] = data.iloc[:, i]
            org = np.unique(data.iloc[:, i]).tolist()
            rep = list(range(len(org)))
            datTraits.replace(to_replace=org, value=rep,
                              inplace=True)
        elif len(np.unique(data.iloc[:, i])) > 2:
            for name in np.unique(data.iloc[:, i]):
                datTraits[name] = data.iloc[:, i]
                org = np.unique(data.iloc[:, i])
                rep = np.repeat(0, len(org))
                rep[np.where(org == name)] = 1
                org = org.tolist()
                rep = rep.tolist()
                datTraits.replace(to_replace=org, value=rep, inplace=True)

    return datTraits


def generate_html_from_df(df: pd.DataFrame, title: str = "Data table", table_id: str = "dataTable", output_path: str = ".",
                          template_path: str = "../flask/app_dev/templates", col_names: List[str] = None) -> str:
    """
    Generate an HTML file from a DataFrame for displaying the analysis table.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the analysis data.
    - title (str): Title for the analysis (used in the HTML title).
    - table_id (str): ID for the DataTable.
    - output_path (str): Path where the HTML file will be saved.
    - template_path (str): Path to the Jinja2 template.
    - col_names (List[str]): List of new column names. If None, keeps the original names.

    Returns:
    - str: Path to the generated HTML file.
    """

    if col_names:
        df.columns = col_names

    df = df.reset_index(drop=True)
    table_html = df.to_html(classes='display dataTable no-border',
                            index=False, border=0, header=True, table_id=table_id)

    # Load the Jinja2 environment
    env = Environment(loader=FileSystemLoader(searchpath=template_path))

    # Load the template
    template = env.get_template("table.html")

    # Render the template with the required dynamic data
    rendered_html = template.render(
        title=title,
        table_html=table_html,
        table_id=table_id
    )

    # Output file path
    output_file = os.path.join(output_path, f"{table_id}.html")

    # Write the rendered HTML to a file
    with open(output_file, "w") as f:
        f.write(rendered_html)

    return output_file


def export_co_expression_network_to_cytoscape(
    tom: Union[pd.DataFrame, List[pd.DataFrame]],
    adata: Union[AnnData, List[AnnData]],
    selected_module_colors: List[str] = None,
    threshold: float = 0.5,
    neighbor_info: Union[set, List[set]] = None
) -> dict:
    """
    Exports the co-expression network(s) in a format compatible with Cytoscape.js.
    
    Parameters:
    - tom (pd.DataFrame or List[pd.DataFrame]): The TOM matrix or list of TOM matrices.
    - adata (AnnData or List[AnnData]): Annotated data object(s).
    - selected_module_colors (List[str]): Module colors to display. If None, all modules are shown.
    - threshold (float): TOM connection threshold for edges.
    - neighbor_info (set or List[set]): Neighbor transcripts for the corresponding TOM(s).
    
    Returns:
    - dict: Dictionary with nodes, edges, and a list of neighbors per TOM.
    """

    if isinstance(tom, list) and isinstance(adata, list):
        if len(tom) != len(adata):
            raise ValueError("The number of TOMs and AnnData objects must be the same.")
        tom_adata_pairs = list(zip(tom, adata))
        neighbor_list = neighbor_info if isinstance(neighbor_info, list) else [neighbor_info] * len(tom)
    elif isinstance(tom, pd.DataFrame) and isinstance(adata, AnnData):
        tom_adata_pairs = [(tom, adata)]
        neighbor_list = [neighbor_info] if neighbor_info is not None else [set()]
    else:
        raise TypeError("Both tom and adata must be either lists or single objects of their respective types.")

    nodes = []
    edges = []
    neighbors_output = []
    organism_id = 0

    for idx, (current_tom, current_adata) in enumerate(tom_adata_pairs):
        gene_metadata = current_adata.var
        if selected_module_colors:
            mod_genes = gene_metadata[gene_metadata['moduleColors'].isin(selected_module_colors)].index
        else:
            mod_genes = gene_metadata.index
        mod_genes = [gene for gene in mod_genes if gene in current_tom.index]
        if 'species' in current_adata.uns:
            organism = current_adata.uns['species']
        else:
            organism = f"Organism_{organism_id}"
            organism_id += 1

        current_neighbors = neighbor_list[idx] if neighbor_list[idx] is not None else set()
        neighbors_output.append({'organism': organism, 'neighbors': list(current_neighbors)})

        for gene in mod_genes:
            node_data = {
                'data': {
                    'id': f"{organism}_{gene}",
                    'gene': gene,
                    'moduleColor': gene_metadata.loc[gene, 'moduleColors'],
                    'ortho_ID': gene_metadata.loc[gene, 'ortho_ID'] if "ortho_ID" in gene_metadata.columns else "",
                    'organism': organism,
                    'name': current_adata.uns['name'],
                    'go_terms': gene_metadata.loc[gene, 'go_terms'] if "go_terms" in gene_metadata.columns else "",
                    'ipr_id': gene_metadata.loc[gene, 'ipr_id'] if "ipr_id" in gene_metadata.columns else ""
                }
            }
            if gene in current_neighbors:
                node_data['data']['is_neighbor'] = True
            node_data['data'] = {k: v if not pd.isna(v) else "" for k, v in node_data['data'].items()}
            nodes.append(node_data)

        for i in range(len(mod_genes)):
            for j in range(i + 1, len(mod_genes)):
                gene_i = mod_genes[i]
                gene_j = mod_genes[j]
                if current_tom.loc[gene_i, gene_j] > threshold:
                    edge_data = {
                        'data': {
                            'source': f"{organism}_{gene_i}",
                            'target': f"{organism}_{gene_j}",
                            'weight': current_tom.loc[gene_i, gene_j],
                            'organism': organism
                        }
                    }
                    edges.append(edge_data)

    return {'nodes': nodes, 'edges': edges, 'neighbors': neighbors_output}


def identify_network_clusters_from_json(network_data: dict, cluster_name: str,
                                        node_threshold: int = None, node_threshold_percent: float = 0.02,
                                        print_info: bool = False) -> Dict[str, str]:
    """
    Assigns nodes to clusters based on the network data in Cytoscape.js format and returns a mapping of nodes to cluster assignments.

    Parameters:
    - network_data (dict): The network data in Cytoscape.js format (contains 'nodes' and 'edges').
    - cluster_name (str): Prefix for cluster assignment.
    - node_threshold (int): Minimum number of nodes required in a cluster to be assigned.
    - node_threshold_percent (float): Percentage of nodes required in a cluster to be assigned.
    - print_info (bool): If True, prints cluster information.

    Returns:
    - Dict[str, str]: Mapping of transcript IDs to cluster assignments.
    """

    if not isinstance(cluster_name, str):
        raise TypeError("cluster_name must be a string")

    # Create an adjacency list from the edges
    adjacency_list = defaultdict(set)
    for edge in network_data['edges']:
        source = edge['data']['source']
        target = edge['data']['target']
        adjacency_list[source].add(target)
        adjacency_list[target].add(source)

    # Find connected components (clusters)
    visited = set()
    clusters = []

    def dfs(node, cluster):
        """Depth-First Search to find all nodes in the same connected component."""
        visited.add(node)
        cluster.append(node)
        for neighbor in adjacency_list[node]:
            if neighbor not in visited:
                dfs(neighbor, cluster)

    for node_data in network_data['nodes']:
        node = node_data['data']['id']
        if node not in visited:
            name = node_data['data']['name']
            cluster = []
            dfs(node, cluster)
            clusters.append((cluster, name))

    # Apply node_threshold based on percentage
    if not node_threshold:
        if not node_threshold_percent:
            raise ValueError(
                "Either node_threshold or node_threshold_percent must be provided")
        node_threshold = int(
            len(network_data['nodes']) * node_threshold_percent)

    # Filter clusters based on the node_threshold
    filtered_clusters = [cluster for cluster in clusters if len(
        cluster[0]) >= node_threshold]

    # Assign clusters to transcripts
    cluster_map = {}
    for new_cluster_id, cluster in enumerate(filtered_clusters, start=1):
        # Get a node name (node_data['data']['name']) from the current cluster
        for transcript in cluster[0]:
            cluster_map[transcript] = f"{cluster[1]}: {cluster_name} {new_cluster_id}"

    # Assign 'No clusters' to transcripts that are not in any cluster
    all_transcripts = {node['data']['id'] for node in network_data['nodes']}
    for transcript in all_transcripts:
        if transcript not in cluster_map:
            cluster_map[transcript] = "No clusters"

    # Optional printing of cluster information
    if print_info:
        cluster_counts = pd.Series(cluster_map).value_counts().sort_index()
        for cluster_id, count in cluster_counts.items():
            cluster_transcripts = [
                t for t, c in cluster_map.items() if c == cluster_id]
            print(
                f"{cluster_id}: {count} transcripts - {', '.join(cluster_transcripts)}")

    return cluster_map


def calculate_tom_metrics(tom: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates TOM-based metrics: Degree, Clustering Coefficient, Connected Components.

    Parameters:
    - tom (pd.DataFrame): The TOM matrix.

    Returns:
    - pd.DataFrame: DataFrame containing TOM-based metrics.
    """

    degrees = (tom > 0).sum(axis=1)
    clustering_coefficient = []

    for gene in tom.index:
        neighbors = tom.loc[gene][tom.loc[gene] > 0].index.tolist()
        if len(neighbors) <= 1:
            clustering_coefficient.append(0)
        else:
            sub_tom = tom.loc[neighbors, neighbors]
            actual_edges = (sub_tom > 0).sum().sum() / 2
            possible_edges = len(neighbors) * (len(neighbors) - 1) / 2
            clustering_coefficient.append(
                actual_edges / possible_edges if possible_edges > 0 else 0)

    connected_components = (degrees > 0).astype(int)

    return pd.DataFrame({
        'degree': degrees,
        'clustering_coefficient': clustering_coefficient,
        'connected_components': connected_components
    })


def calculate_graph_metrics_from_clusters(tom: pd.DataFrame, cluster_info: dict) -> pd.DataFrame:
    """
    Calculates graph-based metrics: Betweenness Centrality, Closeness Centrality, Eccentricity using igraph,
    and ensures all nodes in cluster_info are added to the graph.

    Parameters:
    - tom (pd.DataFrame): The TOM matrix.
    - cluster_info (dict): Dictionary containing cluster assignments.

    Returns:
    - pd.DataFrame: DataFrame containing graph-based metrics.
    """

    # Initialize the iGraph graph with node names
    G = ig.Graph(directed=False)
    G.add_vertices(tom.index.tolist())

    # Add edges with weights (no self-loops)
    edges = [(i, j) for i in range(tom.shape[0])
             for j in range(i + 1, tom.shape[1]) if tom.iat[i, j] > 0]
    weights = [tom.iat[i, j] for i, j in edges]
    G.add_edges(edges)
    G.es['weight'] = weights

    # Check if all nodes in `cluster_info` exist in the graph `G`
    missing_nodes = [node for node in cluster_info if node not in G.vs['name']]
    if missing_nodes:
        # Print missing nodes and filter them out from `cluster_info`
        print(
            f"The following nodes are missing from the graph and will be excluded: {missing_nodes}")
        cluster_info = {node: cluster for node,
                        cluster in cluster_info.items() if node not in missing_nodes}

    # Divide cluster information into subgraphs
    cluster_to_nodes = {}
    for node, cluster in cluster_info.items():
        cluster_to_nodes.setdefault(cluster, []).append(node)

    metrics_list = []
    for cluster, nodes in cluster_to_nodes.items():
        # Create subgraph and ensure it contains nodes
        subgraph = G.subgraph(nodes)
        if not subgraph.vcount():
            continue

        # Calculate metrics for each connected component within the cluster
        components = subgraph.components()
        for component in components:
            component_subgraph = subgraph.subgraph(component)
            num_nodes = len(component_subgraph.vs)

            # Betweenness Centrality, normalized to [0, 1] for the component
            betweenness_centrality = component_subgraph.betweenness(
                weights='weight')
            if num_nodes > 2:
                betweenness_centrality = [
                    b / ((num_nodes - 1) * (num_nodes - 2) / 2) for b in betweenness_centrality]

            # Closeness Centrality and Eccentricity
            closeness_centrality = component_subgraph.closeness()
            eccentricity = [
                int(e) if e is not None else None for e in component_subgraph.eccentricity()]

            # Assign metrics to node names
            for i, node in enumerate(component_subgraph.vs['name']):
                metrics_list.append({
                    'node': node,
                    'cluster': cluster,
                    'betweenness_centrality': betweenness_centrality[i],
                    'closeness_centrality': closeness_centrality[i],
                    'eccentricity': eccentricity[i],
                })

    return pd.DataFrame(metrics_list)


def calculate_all_tom_metrics(tom: Union[pd.DataFrame, List[pd.DataFrame]], adata: Union[AnnData, List[AnnData]],
                              cluster_info: Dict[str, str], tool: str, progress_callback: Callable[[str], None] = None) -> pd.DataFrame:
    """
    Combines TOM-based metrics, cluster-based graph metrics, connectivity, and average TOM value into a single table.
    Can handle a single TOM matrix or a list of TOM matrices.

    Parameters:
    - tom (pd.DataFrame or List[pd.DataFrame]): The TOM matrix or a list of TOM matrices.
    - adata (Union[AnnData, List[AnnData]]): AnnData object(s) to analyze.
    - cluster_info (Dict[str, str]): Dictionary containing cluster assignments for all nodes across all TOMs.
    - tool (str): The tool used to generate the cluster assignments (e.g., "cytoscape").
    - progress_callback (Callable[[str], None]): Callback function to report progress.

    Returns:
    - pd.DataFrame: A combined table containing TOM metrics.
      Each node appears only once, and metrics are calculated based on the TOM it belongs to.
    """

    def calculate_metrics_single(tom_single: pd.DataFrame, adata: AnnData, cluster_info_single: Dict[str, str]) -> pd.DataFrame:
        if tool == "cytoscape":
            # Keep everything after the first "_"
            cluster_info_single = {"_".join(
                key.split("_")[1:]): value for key, value in cluster_info_single.items()}

        if progress_callback:
            progress_callback(
                f"Calculating TOM metrics for dataset {adata.uns['name']}")

        # Fill NaN values with 0
        tom_filled = tom_single.fillna(0)

        # Get species name
        tom_species = adata.uns['species']

        # Calculate TOM-based metrics
        tom_metrics = calculate_tom_metrics(tom_filled)

        # Calculate graph-based metrics from clusters
        graph_metrics = calculate_graph_metrics_from_clusters(
            tom_filled, cluster_info_single)

        # Calculate Connectivity (sum over rows minus self-loops)
        connectivity = tom_filled.sum(axis=1) - tom_filled.values.diagonal()

        # Calculate average TOM value for each node (excluding self-loops)
        average_tom = (tom_filled.sum(
            axis=1) - tom_filled.values.diagonal()) / ((tom_filled > 0).sum(axis=1))
        # Handle NaN values for isolated nodes
        average_tom = average_tom.fillna(0)

        # Set indices for merging
        graph_metrics.set_index('node', inplace=True)
        tom_metrics.index.name = 'node'
        tom_metrics['connectivity'] = connectivity
        tom_metrics['average_tom'] = average_tom

        # Merge TOM metrics and graph metrics
        combined_metrics = pd.merge(
            tom_metrics, graph_metrics, left_index=True, right_index=True)
        combined_metrics['species'] = tom_species

        ortho_ID = adata.var['ortho_ID'] if 'ortho_ID' in adata.var.columns else pd.Series(
            [""] * len(adata.var), index=adata.var.index)
        combined_metrics['ortho_ID'] = combined_metrics.index.map(ortho_ID)

        # Set ortho_ID as the second column
        combined_metrics = combined_metrics[[
            'ortho_ID'] + [col for col in combined_metrics.columns if col != 'ortho_ID']]
        combined_metrics['ortho_ID'] = combined_metrics['ortho_ID'].cat.rename_categories({
                                                                                          '': 'No orthogroup'})
        # Add the new category 'No orthogroup' to the Categorical column
        combined_metrics['ortho_ID'] = combined_metrics['ortho_ID'].cat.add_categories('No orthogroup')

        # Replace NaN values with 'No orthogroup'
        combined_metrics['ortho_ID'] = combined_metrics['ortho_ID'].fillna('No orthogroup')


        # Set Species as the first column
        combined_metrics = combined_metrics[[
            'species'] + [col for col in combined_metrics.columns if col != 'species']]

        # Set node as index
        combined_metrics.index.name = 'Node'

        # Rename columns for consistency
        combined_metrics.columns = [
            'Species', 'Orthogroup', 'Degree', 'Clustering Coefficient', 'Connected Components', 'Connectivity',
            'Average TOM Value', 'Cluster', 'Betweenness Centrality', 'Closeness Centrality', 'Eccentricity'
        ]

        # Drop columns Connectivity, Average TOM Value, and Connected Components
        combined_metrics.drop(
            columns=['Connectivity', 'Average TOM Value', 'Connected Components'], inplace=True)

        return combined_metrics

    if isinstance(tom, list):
        if not isinstance(cluster_info, dict):
            raise ValueError(
                "When 'tom' is a list, 'cluster_info' must be a single dictionary containing cluster assignments for all nodes.")
        if not isinstance(adata, list):
            raise ValueError(
                "When 'tom' is a list, 'adata' must be a list of AnnData objects.")
        if not len(tom) == len(adata):
            raise ValueError(
                "The number of TOMs and AnnData objects must be the same.")

        all_metrics = []

        for i, tom_single in enumerate(tom):
            # Check if tom_single is empty
            if tom_single.empty:
                continue

            # Get nodes in this TOM
            nodes = tom_single.index

            # Add species prefix to nodes
            nodes = [f"{adata[i].uns['species']}_{node}" for node in nodes]

            # Extract cluster_info for these nodes
            cluster_info_single = {
                node: cluster_info[node] for node in nodes if node in cluster_info}

            # Calculate metrics for the current TOM
            metrics_single = calculate_metrics_single(
                tom_single, adata[i], cluster_info_single)

            # Append to the list
            all_metrics.append(metrics_single)

        # Concatenate all metrics into a single DataFrame
        combined_all_metrics = pd.concat(all_metrics)

        return combined_all_metrics

    elif isinstance(tom, pd.DataFrame):
        if not isinstance(cluster_info, dict):
            raise ValueError(
                "When 'tom' is a single DataFrame, 'cluster_info' must be a single dictionary containing cluster assignments.")

        # Calculate metrics for the single TOM
        combined_metrics = calculate_metrics_single(tom, adata, cluster_info)

        return combined_metrics

    else:
        raise TypeError(
            "'tom' must be either a pandas DataFrame or a list of DataFrames.")


def add_go_terms_to_adata(adata: AnnData, go_mapping_file: str) -> AnnData:
    """
    Add a 'go_terms' column to adata.var based on a mapping file.

    Parameters:
    - adata (AnnData): AnnData object with var index matching transcript IDs.
    - go_mapping_file (str): Path to the file containing GO term mappings.

    Returns:
    - AnnData: Updated AnnData object with a new 'go_terms' column in var.
    """

    if not os.path.exists(go_mapping_file):
        print(
            f"Warning: GO mapping file '{go_mapping_file}' not found. Filling 'go_terms' with empty strings.")
        adata.var["go_terms"] = ''
        return adata

    # Load GO term mapping file into a DataFrame
    go_df = pd.read_csv(go_mapping_file, sep="\t",
                        header=None, names=["transcript", "go_terms"])

    # Ensure GO terms are separated by commas instead of spaces
    go_df["go_terms"] = go_df["go_terms"].str.replace(r"\s+", ",", regex=True)

    # Group GO terms by transcript and join them as comma-separated strings
    go_df = go_df.groupby("transcript")["go_terms"].apply(
        lambda x: ",".join(x)).reset_index()

    # Create a dictionary for fast lookup
    go_dict = go_df.set_index("transcript")["go_terms"].to_dict()

    # Map GO terms to adata.var.index
    adata.var["go_terms"] = adata.var.index.map(go_dict)

    # Replace missing mappings with empty strings
    adata.var["go_terms"] = adata.var["go_terms"].fillna('')

    return adata


def add_ipr_columns(adata: AnnData, ipr_mapping_file: str) -> AnnData:
    """
    Add 'ipr_id' and 'ipr_desc' columns to adata.var based on an IPR mapping file.

    Parameters:
    - adata (AnnData): AnnData object with var index matching transcript IDs.
    - ipr_mapping_file (str): Path to the file containing IPR ID and description mappings.

    Returns:
    - AnnData: Updated AnnData object with new 'ipr_id' and 'ipr_desc' columns.
    """

    if not os.path.exists(ipr_mapping_file):
        print(
            f"Warning: IPR mapping file '{ipr_mapping_file}' not found. Filling 'ipr_id' and 'ipr_desc' with empty strings.")
        adata.var["ipr_id"] = ''
        adata.var["ipr_desc"] = ''
        return adata

    # Load the IPR mapping file
    ipr_df = pd.read_csv(ipr_mapping_file, sep="\t", header=None, names=[
                         "transcript", "ipr_id", "ipr_desc"])

    # Group IPR IDs and descriptions by transcript
    grouped_ipr = ipr_df.groupby("transcript").agg({
        "ipr_id": lambda x: ",".join(x),
        "ipr_desc": lambda x: ",".join(x)
    }).reset_index()

    # Create dictionaries for fast mapping
    ipr_dict_ids = grouped_ipr.set_index("transcript")["ipr_id"].to_dict()
    ipr_dict_desc = grouped_ipr.set_index("transcript")["ipr_desc"].to_dict()

    # Map the IPR data to adata.var
    adata.var["ipr_id"] = adata.var.index.map(ipr_dict_ids)
    adata.var["ipr_desc"] = adata.var.index.map(ipr_dict_desc)

    # Fill missing values with empty strings
    adata.var["ipr_id"] = adata.var["ipr_id"].fillna('')
    adata.var["ipr_desc"] = adata.var["ipr_desc"].fillna('')

    return adata


def add_combined_column(adata: AnnData, columns: List[str] = None, new_column_name: str = "Combined_Trait", drop_others: bool = False) -> None:
    """
    Adds a new column to adata.obs by combining values from specified or all columns.
    Optionally drops all other columns after adding the combined column.

    Parameters:
    - adata (AnnData): The AnnData object containing the obs DataFrame.
    - columns (List[str]): The columns to combine. If None, all columns are used.
    - new_column_name (str): The name of the new combined column.
    - drop_others (bool): Whether to drop all other columns after adding the combined column.
    """

    # Use all columns if none are specified
    if columns is None:
        columns = adata.obs.columns.tolist()
    else:
        # Ensure specified columns exist
        for col in columns:
            if col not in adata.obs.columns:
                raise ValueError(f"Column '{col}' not found in adata.obs")

    # Combine columns into the new column
    adata.obs[new_column_name] = adata.obs[columns].apply(
        lambda row: " ".join(map(str, row)), axis=1
    )

    # Drop all other columns if specified
    if drop_others:
        adata.obs = adata.obs[[new_column_name]]
        print(
            f"Added combined column '{new_column_name}' and dropped all other columns")
    else:
        print(
            f"Added combined column '{new_column_name}', other columns retained")
