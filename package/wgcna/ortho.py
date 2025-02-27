import gzip
import os
import pandas as pd
import numpy as np
from anndata import AnnData
from typing import List, Union, Dict, Tuple
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from concurrent.futures import ProcessPoolExecutor, as_completed

##########
"""
USE FLAT FILES
"""
##########


def get_transcripts_by_hog(hog_file: str, hog_ids: List[str]) -> Dict[str, List[str]]:
    """
    Extract transcript lists for given HOG IDs from the HOG file.

    Parameters:
    - hog_file (str): Path to the HOG file.
    - hog_ids (List[str]): List of specific HOG IDs.

    Returns:
    - Dict[str, List[str]]: A dictionary with plants as keys and lists of transcripts as values.
    """

    transcripts_by_plant = {}
    hog_ids = set(hog_ids)
    with open(hog_file, 'r') as file:
        headers = file.readline().strip().split('\t')
        for line in file:
            parts = line.strip().split('\t')
            if parts[0] in hog_ids:
                if len(parts) > 3:
                    for i, header in enumerate(headers[4:], start=4):
                        if i < len(parts) and parts[i].strip():
                            transcripts = parts[i].strip().split(', ')
                            if header in transcripts_by_plant:
                                transcripts_by_plant[header].update(
                                    transcripts)
                            else:
                                transcripts_by_plant[header] = set(transcripts)
                        else:
                            transcripts_by_plant.setdefault(header, set())

    return {k: list(v) for k, v in transcripts_by_plant.items()}


def create_transcript_df(transcripts: List[str], gaf_file_path: str) -> pd.DataFrame:
    """
    Create a DataFrame for transcripts with their GO Terms and Gene Names.

    Parameters:
    - transcripts (List[str]): List of transcript IDs.
    - gaf_file_path (str): Path to the GAF.gz file.

    Returns:
    - pd.DataFrame: DataFrame with transcripts as the index, and GO Terms and Gene Names as columns.
    """

    result = {transcript: {'GO_Terms': [], 'Gene_Name': None}
              for transcript in transcripts}

    with gzip.open(gaf_file_path, 'rt') as file:
        for line in file:
            if line.startswith('!'):
                continue  # Skip comment lines
            parts = line.strip().split('\t')
            gene_id = parts[1]  # Extract the gene ID
            go_term = parts[4]  # Extract the GO Term
            gene_name = parts[2]  # Extract the Gene Name

            if gene_id in transcripts:
                if result[gene_id]['Gene_Name'] is None:
                    result[gene_id]['Gene_Name'] = gene_name

                result[gene_id]['GO_Terms'].append(go_term)

    for transcript in result:
        result[transcript]['GO_Terms'] = ','.join(
            set(result[transcript]['GO_Terms']))

    df = pd.DataFrame.from_dict(result, orient='index')

    return df


def create_gaf_dict(directory: str, suffix: str = ".annotation.goslim_plant.gaf.gz") -> Dict[str, str]:
    """
    Create a dictionary for GAF file paths based on the files in a given directory.

    Parameters:
    - directory (str): Path to the directory containing the GAF files.
    - suffix (str): Suffix for GAF files. Defaults to ".annotation.goslim_plant.gaf.gz".

    Returns:
    - Dict[str, str]: A dictionary with plant names as keys and paths to GAF files as values.
    """

    gaf_files_dict = {}

    for filename in os.listdir(directory):
        if filename.endswith('.gaf.gz') and suffix in filename:
            plant_name = filename.split('.')[0]
            gaf_files_dict[plant_name] = os.path.join(directory, filename)

    return gaf_files_dict


def process_hog_to_dfs(hog_ids: List[str], hog_file: str, gaf_files_by_plant: Dict[str, str], verbose: bool = True) -> Dict[str, pd.DataFrame]:
    """
    Process a list of HOG IDs to create dataframes with transcripts, GO terms, and gene names for each plant.

    Parameters:
    - hog_ids (List[str]): List of HOG IDs to process.
    - hog_file (str): Path to the HOG file.
    - gaf_files_by_plant (Dict[str, str]): Dictionary mapping plant names to paths of GAF files.
    - verbose (bool): Flag to control the output of print statements for processing information.

    Returns:
    - Dict[str, pd.DataFrame]: A dictionary with plant names as keys and corresponding dataframes as values.
    """

    all_go_terms = set()
    all_gene_names = set()
    dataframes = {}

    for hog_id in hog_ids:
        transcripts_by_plant = get_transcripts_by_hog(hog_file, [hog_id])
        for plant, transcripts in transcripts_by_plant.items():
            if transcripts:
                if verbose:
                    print(
                        f"Processing {plant} with {len(transcripts)} transcripts for HOG ID {hog_id}.")
                df = create_transcript_df(
                    transcripts, gaf_files_by_plant[plant])
                if plant not in dataframes:
                    dataframes[plant] = df
                else:
                    dataframes[plant] = pd.concat([dataframes[plant], df])
                all_go_terms.update(','.join(df['GO_Terms']).split(','))
                all_gene_names.update(df['Gene_Name'].dropna().unique())

    if verbose:
        print(
            f"Total unique GO Terms: {len(all_go_terms)} - {list(all_go_terms)}")
        print(
            f"Total unique Gene Names: {len(all_gene_names)} - {list(all_gene_names)}")

    return dataframes


def find_hogs_in_n_species(hog_file: str, n: int) -> List[str]:
    """
    Returns a list of HOG IDs that occur in exactly n different species.

    Parameters:
    - hog_file (str): Path to the HOG file.
    - n (int): The exact number of species in which the HOG must appear.

    Returns:
    - list: List of HOG IDs that occur in exactly n species.
    """

    hog_counts = {}

    with open(hog_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            hog_id = parts[0]
            species_count = 0

            for i in range(3, len(parts)):  # Skip the first four columns
                if parts[i].strip():  # Check if the column is not empty
                    species_count += 1

            if species_count in hog_counts:
                hog_counts[species_count].append(hog_id)
            else:
                hog_counts[species_count] = [hog_id]

    # Returns an empty list if no HOGs meet the criteria
    return hog_counts.get(n, [])


##########
"""
USE ADATAS
"""
##########


def add_ortho_count(adata_or_list: Union[AnnData, List[AnnData]], ortho_id_col: str = 'ortho_id', ortho_count_col: str = 'ortho_count') -> None:
    """
    Adds and calculates the 'ortho_count' column to the provided AnnData object(s),
    if it does not already exist. Empty strings in 'ortho_id' are not counted.

    Parameters:
    - adata_or_list (Union[AnnData, List[AnnData]]): An AnnData object or a list of AnnData objects.
    - ortho_id_col (str): The column name for ortholog IDs. 
    - ortho_count_col (str): The column name for ortholog counts. 
    """

    def process_adata(adata):
        if ortho_count_col in adata.var:
            adata.var.drop(columns=ortho_count_col, errors='ignore', inplace=True)

        # Filter out empty strings before counting
        valid_ortho_ids = adata.var[ortho_id_col][adata.var[ortho_id_col] != '']
        ortho_counts = valid_ortho_ids.value_counts()
        adata.var[ortho_count_col] = adata.var[ortho_id_col].map(ortho_counts).fillna(0).astype(int)

    if isinstance(adata_or_list, list):
        for adata in adata_or_list:
            process_adata(adata)
    else:
        process_adata(adata_or_list)


def update_ortho_ID(adatas: Union[AnnData, List[AnnData]], node: int = 0, ortho_id_col: str = 'ortho_id', node_prefix: str = 'N', node_suffix: str = '_ortho') -> None:
    """
    Updates the 'ortho_id' column in the adata.var dataframe to the values corresponding to the specified node.

    Parameters:
    - adatas (Union[AnnData, List[AnnData]]): An AnnData object or a list of AnnData objects.
    - node (int): The node level to select the corresponding ortho_ids. 
    - ortho_id_col (str): The column name for ortholog IDs. 
    - node_prefix (str): The prefix for node columns. 
    - node_suffix (str): The suffix for node columns. 
    """

    def set_ortho_id(adata: AnnData, node: int, ortho_id_col: str, node_prefix: str, node_suffix: str):
        node_key = f"{node_prefix}{node}{node_suffix}"
        if node_key in adata.var.columns:
            adata.var[ortho_id_col] = adata.var[node_key]
        else:
            raise ValueError(
                f"Node {node} is not available in the adata.var dataframe.")

    if isinstance(adatas, AnnData):
        set_ortho_id(adatas, node, ortho_id_col, node_prefix, node_suffix)
        add_ortho_count(adatas, ortho_id_col=ortho_id_col)
    else:
        for adata in adatas:
            set_ortho_id(adata, node, ortho_id_col, node_prefix, node_suffix)
            add_ortho_count(adata, ortho_id_col=ortho_id_col)


def get_ortho_df(adatas: List[AnnData], 
                 transcripts: Dict[str, List[str]] = None, 
                 ortho_id_col: str = 'ortho_id', 
                 ortho_count_col: str = 'ortho_count', 
                 count_col: str = 'Count') -> pd.DataFrame:
    """
    Calculates a DataFrame with ortho_id as the index and ortho_count for each species,
    including a Count column summarizing the ortho_count across species. Optionally filters
    the data based on a dictionary of transcripts by species.

    Parameters:
    - adatas (List[AnnData]): A list of AnnData objects.
    - transcripts (Dict[str, List[str]]): Dictionary with species as keys and lists of transcript names to filter the data.
    - ortho_id_col (str): Column name for ortholog IDs. 
    - ortho_count_col (str): Column name for ortholog counts. 
    - count_col (str): Column name for the total count. 

    Returns:
    - pd.DataFrame: A DataFrame containing ortho_id and ortho_count for each species,
      and a Count column.
    """

    species_names = [adata.uns['species'] for adata in adatas]
    orthogroups = set()

    for adata in adatas:
        species = adata.uns['species']

        # Check if species-specific transcripts are provided and filter accordingly
        if transcripts is not None and species in transcripts:
            # Filter adata.var to only include specified transcripts for the species
            adata = adata[:, adata.var_names.isin(transcripts[species])]

        # Collect unique orthogroups across all datasets
        orthogroups.update(adata.var[ortho_id_col])

    orthogroups = list(orthogroups)

    # Create an empty DataFrame with Orthogroups as the index
    df = pd.DataFrame(index=orthogroups, columns=species_names + [count_col])

    # Fill the DataFrame with ortho counts per species
    for adata in adatas:
        species = adata.uns['species']

        # Filter adata.var again to only include specified transcripts for this species, if provided
        if transcripts is not None and species in transcripts:
            adata = adata[:, adata.var_names.isin(transcripts[species])]

        # Extract unique orthogroups and their counts
        orthogroup_counts = adata.var.drop_duplicates(subset=[ortho_id_col])
        for _, row in orthogroup_counts.iterrows():
            orthogroup = row[ortho_id_col]
            count = row[ortho_count_col]
            df.at[orthogroup, species] = count

    # Remove any rows where the index is NaN
    df = df[df.index.notna()]
    df.fillna(0, inplace=True)

    # Calculate the total ortho_count across all species and add it to the Count column
    df[count_col] = df.iloc[:, :-1].sum(axis=1)
    df.sort_values(by=count_col, ascending=False, inplace=True)

    return df


def get_filtered_ortho_df(filter_df: pd.DataFrame, 
                          transcripts: List[str] = None, 
                          ortho_id_col: str = 'ortho_id', 
                          species_col: str = 'species', 
                          transcript_col: str = 'Transcript', 
                          count_col: str = 'Count') -> pd.DataFrame:
    """
    Creates a DataFrame with ortho_id as the index and species as columns, with each cell containing
    the count of transcripts per species for each ortho_id, and includes a Count column summarizing 
    the total counts across species.

    Parameters:
    - filter_df (pd.DataFrame): A DataFrame with a MultiIndex (species, Transcript) and containing ortho_id.
    - transcripts (List[str]): A list of transcript names to filter the DataFrame.
    - ortho_id_col (str): Column name for ortholog IDs. 
    - species_col (str): Column name for species. 
    - transcript_col (str): Column name for transcripts. 
    - count_col (str): Column name for the total count. 

    Returns:
    - pd.DataFrame: A DataFrame containing ortho_id as the index, species as columns, 
      and a Count column summarizing the transcript counts across species.
    """

    if transcripts is not None:
        filter_df = filter_df[filter_df.index.get_level_values(transcript_col).isin(transcripts)]

    filter_df = filter_df[filter_df[ortho_id_col].notna() & (filter_df[ortho_id_col] != "")]
    grouped = filter_df.groupby([ortho_id_col, species_col]).size().unstack(fill_value=0)
    grouped[count_col] = grouped.sum(axis=1)
    grouped.sort_values(by=count_col, ascending=False, inplace=True)
    grouped = grouped.astype(int)
    grouped.index.name = None
    grouped.columns.name = None

    return grouped


def extract_go_terms(df: pd.DataFrame) -> List[str]:
    """
    Extract comma-separated GO terms from the 'go_terms' column
    and return a list of strings without NaN values.

    - df (pd.DataFrame): DataFrame containing the 'go_terms' column

    Returns:
    - List[str]: List of GO terms.
    """

    go_terms = df['go_terms'].dropna().str.split(',')

    return list(set([term for sublist in go_terms for term in sublist]))


def transcript_ortho_browser(
    query: Union[str, int, List[str], List[int], None],
    adatas: Union[AnnData, List[AnnData]], 
    possible_columns: List[str] = ['module_colors', 'module_labels', 'go_terms', 'ortho_id'],
    species_key: str = 'species',
    module_colors_col: str = 'module_colors',
    ortho_id_col: str = 'ortho_id',
    module_count_col: str = 'module_count',
    ortho_count_col: str = 'ortho_count',
    unique_modules_col: str = 'unique_modules',
    unique_go_terms_col: str = 'unique_go_terms',
    total_counts_col: str = 'total_counts'
) -> pd.DataFrame:
    """
    Implements a Transcript/Ortho-Browser. Given an ortho ID, GO term, transcript, or moduleColor (or a list of them),
    this function returns a DataFrame with the relevant information.

    Parameters:
    - query (Union[str, int, List[str], List[int], None]): The ortho ID, GO term, transcript, moduleColor, or moduleLabel 
      (or a list of them) to query. If None or empty string, no filtering is applied.
    - adatas (Union[AnnData, List[AnnData]]): A list of AnnData objects or a single AnnData object.
    - possible_columns (List[str]): List of possible columns to include in the result DataFrame.
    - species_key (str): Key for species information in AnnData.uns.
    - module_colors_col (str): Column name for module colors.
    - ortho_id_col (str): Column name for ortholog IDs.
    - module_count_col (str): Column name for module counts.
    - ortho_count_col (str): Column name for ortholog counts.
    - unique_modules_col (str): Column name for unique modules.
    - unique_go_terms_col (str): Column name for unique GO terms.
    - total_counts_col (str): Column name for total counts.

    Returns:
    - pd.DataFrame: A DataFrame with the relevant information for the query.
    """

    # Ensure adatas is a list
    if isinstance(adatas, AnnData):
        adatas = [adatas]

    # Normalize the query to a list
    if isinstance(query, (str, int)):
        query = [query]

    results = []

    for adata in adatas:
        species = adata.uns.get(species_key, 'unknown')

        # Select available columns in the current adata.var
        columns = [col for col in possible_columns if col in adata.var.columns]

        # Handle case where no query is provided (return all rows)
        if query is None or (isinstance(query, list) and all(q == "" for q in query)):
            df = adata.var[columns].copy()
        else:
            # Filter by query
            df_list = []
            for q in query:
                q_str = str(q)
                if q in adata.var.index:
                    df_list.append(adata.var.loc[[q], columns])
                elif any(col in adata.var.columns for col in possible_columns):
                    for col in possible_columns:
                        if col in adata.var.columns and any(adata.var[col].astype(str).str.contains(q_str, na=False)):
                            df_list.append(adata.var[adata.var[col].astype(
                                str).str.contains(q_str, na=False)][columns])

            # Combine filtered DataFrames
            df = pd.concat(df_list, axis=0) if df_list else pd.DataFrame(
                columns=columns)

        if not df.empty:
            df[species_key] = species

            # Add 'module_count' if 'module_colors' exists
            if module_colors_col in df.columns:
                module_counts = df[module_colors_col].value_counts().to_dict()
                df[module_count_col] = df[module_colors_col].map(module_counts)

            # Add 'ortho_count' if 'ortho_id' exists
            if ortho_id_col in df.columns:
                valid_ortho_ids = adata.var[ortho_id_col][adata.var[ortho_id_col] != '']
                ortho_counts = valid_ortho_ids.value_counts()
                df[ortho_count_col] = df[ortho_id_col].map(
                    ortho_counts).fillna(0).astype(int)

            results.append(df)

    if results:
        # Combine all results
        result_df = pd.concat(results)
        result_df.index.name = "transcript"

        # Set species and transcript as multi-index
        result_df.set_index([species_key, result_df.index], inplace=True)

        # Calculate unique modules and GO terms per species (if applicable)
        if module_colors_col in result_df.columns:
            unique_modules_per_species = result_df.groupby(
                species_key)[module_colors_col].nunique().to_dict()
            result_df[unique_modules_col] = result_df.index.get_level_values(
                species_key).map(unique_modules_per_species)

        if 'go_terms' in result_df.columns:
            unique_go_terms_per_species = result_df.groupby(
                species_key)['go_terms'].nunique().to_dict()
            result_df[unique_go_terms_col] = result_df.index.get_level_values(
                species_key).map(unique_go_terms_per_species)

        # Ensure numeric types for 'total_counts' if it exists
        if total_counts_col in result_df.columns:
            result_df[total_counts_col] = result_df[total_counts_col].fillna(
                0).astype(int)

        # Standardize column names to lowercase
        result_df.columns = [col.lower() for col in result_df.columns]
        result_df.index.names = [col.lower() for col in result_df.index.names]

        return result_df
    else:
        # Return an empty DataFrame with all possible columns
        empty_columns = [col.lower() for col in possible_columns + [species_key,
                                                                    module_count_col, ortho_count_col, unique_modules_col, unique_go_terms_col]]

        return pd.DataFrame(columns=empty_columns)


def prepare_dataframe(data: List[Dict[str, Dict[str, str]]], 
                      cluster_info: Dict[str, str] = None, 
                      data_key: str = 'data', 
                      id_key: str = 'id', 
                      ortho_id_key: str = 'ortho_id', 
                      organism_key: str = 'species', 
                      cluster_col: str = 'Cluster', 
                      no_clusters_label: str = 'No clusters', 
                      no_orthogroup_label: str = 'No orthogroup') -> pd.DataFrame:
    """
    Prepares a DataFrame from the raw transcript data, mapping cluster information and handling
    empty ortho_id values.

    Parameters:
    - data (List[Dict[str, Dict[str, str]]]): A list of dictionaries containing transcript information.
      Each dictionary should have a 'data' key with values for 'id', 'ortho_id', 'organism', and other metadata.
    - cluster_info (Dict[str, str]): A dictionary mapping node ids to cluster names.
    - data_key (str): Key to access the data in the input dictionaries. 
    - id_key (str): Key for the transcript ID in the data. 
    - ortho_id_key (str): Key for the ortholog ID in the data. 
    - organism_key (str): Key for the organism in the data. 
    - cluster_col (str): Column name for cluster information in the output DataFrame. 
    - no_clusters_label (str): Label for entries with no cluster information. 
    - no_orthogroup_label (str): Label for entries with no orthogroup information. 

    Returns:
    - pd.DataFrame: A processed DataFrame with cluster information mapped and 'No orthogroup'
      labels applied to empty ortho_id values.
    """
    # Convert data to DataFrame
    df = pd.DataFrame([entry[data_key] for entry in data])

    # Map cluster information to the DataFrame; default to 'No clusters' if id not in cluster_info
    if cluster_info:
        df[cluster_col] = df[id_key].map(cluster_info).fillna(no_clusters_label)

    # Replace empty and NaN ortho_id values with 'No orthogroup'
    df[ortho_id_key] = df[ortho_id_key].replace('', no_orthogroup_label).fillna(no_orthogroup_label)

    # Drop all rows where "Cluster" is "No clusters"
    df = df[df[cluster_col] != no_clusters_label]

    return df


def summarize_orthogroups(df: pd.DataFrame, 
                          ortho_id_col: str = 'ortho_id', 
                          id_col: str = 'id', 
                          organism_col: str = 'species', 
                          cluster_col: str = 'Cluster', 
                          transcripts_col: str = 'Transcripts', 
                          species_col: str = 'Species', 
                          clusters_col: str = 'Clusters', 
                          orthogroup_index_name: str = 'Orthogroup') -> pd.DataFrame:
    """
    Summarizes transcript data by orthogroups, calculating the number of transcripts, unique organisms, 
    and unique clusters per orthogroup.

    Parameters:
    - df (pd.DataFrame): The prepared DataFrame with columns for orthogroup ID, transcript ID, organism, and cluster.
    - ortho_id_col (str): Column name for orthogroup IDs. 
    - id_col (str): Column name for transcript IDs. 
    - organism_col (str): Column name for organism names. 
    - cluster_col (str): Column name for cluster names. 
    - transcripts_col (str): Column name for the number of transcripts. 
    - species_col (str): Column name for the number of unique organisms. 
    - clusters_col (str): Column name for the number of unique clusters. 
    - orthogroup_index_name (str): Name for the orthogroup index. 

    Returns:
    - pd.DataFrame: A DataFrame indexed by orthogroup with columns for the number of transcripts, unique organisms, and unique clusters.
    """

    # Group by ortho_id_col and calculate summary metrics
    df_summary = df.groupby(ortho_id_col).agg(
        # Count the number of transcripts (ids)
        **{transcripts_col: (id_col, 'count')},
        # Count the number of unique organisms
        **{species_col: (organism_col, 'nunique')},
        # Count the number of unique clusters
        **{clusters_col: (cluster_col, 'nunique')}
    ).reset_index()

    # Set ortho_id_col as the index and rename it to orthogroup_index_name
    df_summary.set_index(ortho_id_col, inplace=True)
    df_summary.index.name = orthogroup_index_name

    return df_summary


def calculate_jaccard_matrix(df: pd.DataFrame, cluster_col: str = "Cluster", ortho_id_col: str = "ortho_id", organism_col: str = "species") -> pd.DataFrame:
    """
    Calculates a square matrix containing the Jaccard Index between clusters
    based on shared 'ortho_id' values, including species information in the index.

    Parameters:
    - df (pd.DataFrame): DataFrame containing at least the cluster, ortho_id, and organism columns.
    - cluster_col (str): Name of the column representing cluster labels.
    - ortho_id_col (str): Name of the column representing ortholog IDs.
    - organism_col (str): Name of the column representing species information.

    Returns:
    - pd.DataFrame: A symmetric matrix where each entry [i, j] represents the Jaccard Index
      between clusters i and j based on their shared 'ortho_id's, indexed by species.
    """
    # Extract unique species-cluster combinations
    df['species_cluster'] = df[organism_col] + \
        ": " + df[cluster_col].astype(str)
    clusters = df['species_cluster'].unique()

    # Create dictionary for ortho_ids per species-cluster
    cluster_orthos = {cluster: set(
        df[df['species_cluster'] == cluster][ortho_id_col]) for cluster in clusters}

    # Initialize an empty matrix for Jaccard indices with species-cluster as index and columns
    jaccard_matrix = pd.DataFrame(
        np.zeros((len(clusters), len(clusters))), index=clusters, columns=clusters)

    # Calculate Jaccard index for each cluster pair
    for i, cluster1 in enumerate(clusters):
        for j, cluster2 in enumerate(clusters):
            if i <= j:  # Calculate only the upper triangular matrix
                intersection = cluster_orthos[cluster1].intersection(
                    cluster_orthos[cluster2])
                union = cluster_orthos[cluster1].union(
                    cluster_orthos[cluster2])
                jaccard_index = len(intersection) / \
                    len(union) if len(union) > 0 else 0
                jaccard_matrix.loc[cluster1, cluster2] = jaccard_index
                # Symmetric entry
                jaccard_matrix.loc[cluster2, cluster1] = jaccard_index

    # Set the diagonal values to 1 for self-similarity
    np.fill_diagonal(jaccard_matrix.values, 1)

    return jaccard_matrix


###############
"""
Compare species
"""
###############

# START COMMON HOGS AND JACCARD KOEFFICIENT


def calculate_module_similarity(adata_list: List[AnnData], 
                                module_col: str = 'module_colors', 
                                ortho_id_col: str = 'ortho_id', 
                                species_key: str = 'species') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compares modules across multiple AnnData objects using HOGs.

    Parameters:
    - adata_list (List[AnnData]): List of AnnData objects representing different plant species.
    - module_col (str): Column name for module colors. 
    - ortho_id_col (str): Column name for ortholog IDs. 
    - species_key (str): Key for species information in AnnData.uns. 

    Returns:
    - Tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames, one containing the number of common HOGs per module 
      and another containing the Jaccard coefficient of HOGs per module.
    """

    module_hogs = {}

    # Extract HOGs for each module in each species
    for adata in adata_list:
        species = adata.uns[species_key]
        for module in adata.var[module_col].unique():
            hogs = set(adata.var[adata.var[module_col] == module][ortho_id_col])
            module_hogs[f"{module}_{species}"] = hogs

    # Initialize matrices
    module_names = list(module_hogs.keys())
    num_modules = len(module_names)

    common_hogs_matrix = pd.DataFrame(
        np.zeros((num_modules, num_modules)), index=module_names, columns=module_names)
    jaccard_matrix = pd.DataFrame(
        np.zeros((num_modules, num_modules)), index=module_names, columns=module_names)

    # Calculate common HOGs and Jaccard coefficient
    for i, mod1 in enumerate(module_names):
        for j, mod2 in enumerate(module_names):
            if i <= j:
                hogs1 = module_hogs[mod1]
                hogs2 = module_hogs[mod2]

                common_hogs = hogs1.intersection(hogs2)
                union_hogs = hogs1.union(hogs2)

                common_hogs_matrix.at[mod1, mod2] = common_hogs_matrix.at[mod2, mod1] = len(common_hogs)
                jaccard_index = len(common_hogs) / len(union_hogs) if union_hogs else 0
                jaccard_matrix.at[mod1, mod2] = jaccard_matrix.at[mod2, mod1] = jaccard_index

    return common_hogs_matrix, jaccard_matrix

# END COMMON HOGS AND JACCARD KOEFFICIENT

# START CORRELATION OF EIGENGENES


def aggregate_eigengenes(adata: AnnData, eigengenes_key: str = 'eigengenes', ortho_id_col: str = 'ortho_id') -> pd.DataFrame:
    """
    Aggregates the eigengenes by ortho_id by taking the mean of eigengenes.

    Parameters:
    - adata (AnnData): The AnnData object.
    - eigengenes_key (str): Key for accessing eigengenes in adata.varm. 
    - ortho_id_col (str): Column name for ortholog IDs. 

    Returns:
    - pd.DataFrame: Aggregated eigengenes DataFrame.
    """

    eigengenes_df = pd.DataFrame(
        adata.varm[eigengenes_key], index=adata.var_names)
    eigengenes_df[ortho_id_col] = adata.var[ortho_id_col].values
    aggregated_df = eigengenes_df.groupby(ortho_id_col).mean()

    return aggregated_df


def calculate_correlation_pair(col1: str, col2: str, combined_df: pd.DataFrame) -> Tuple[Tuple[str, str], float, float]:
    """
    Calculate the Spearman correlation and p-value for a pair of columns.

    Parameters:
    - col1 (str): The first column name.
    - col2 (str): The second column name.
    - combined_df (pd.DataFrame): The combined DataFrame with aggregated eigengenes.

    Returns:
    - Tuple[Tuple[str, str], float, float]: Column pair, correlation, and p-value.
    """

    corr, p_value = spearmanr(combined_df[col1], combined_df[col2])

    return (col1, col2), corr, p_value


def calculate_correlation_matrix(adatas: List[AnnData], 
                                 correction_method: str = 'fdr_bh', 
                                 ortho_id_col: str = 'ortho_id', 
                                 species_key: str = 'species', 
                                 eigengenes_key: str = 'eigengenes', 
                                 kme_prefix: str = 'kME') -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Calculates the correlation matrix of aggregated eigengenes between multiple AnnData objects and returns p-value matrix.

    Parameters:
    - adatas (List[AnnData]): List of AnnData objects.
    - correction_method (str): Method to use for multiple testing correction.
                               Options include 'bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg',
                               'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky'.
    - ortho_id_col (str): Column name for ortholog IDs. 
    - species_key (str): Key for species information in AnnData.uns. 
    - eigengenes_key (str): Key for accessing eigengenes in adata.varm. 
    - kme_prefix (str): Prefix for kME columns. 

    Returns:
    - Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: Correlation matrix, p-value matrix, and combined DataFrame.
    """

    all_hogs = set()
    for adata in adatas:
        all_hogs.update(adata.var[ortho_id_col])

    aggregated_dfs = []
    species_names = []

    for adata in adatas:
        aggregated_df = aggregate_eigengenes(adata, eigengenes_key=eigengenes_key, ortho_id_col=ortho_id_col)
        species = adata.uns[species_key]
        aggregated_df.columns = [f"{col}_{species}" for col in aggregated_df.columns]
        species_names.append(species)
        aggregated_dfs.append(aggregated_df)

    combined_df = pd.concat(aggregated_dfs, axis=1, join='outer').fillna(0)

    correlation_matrix = pd.DataFrame(np.zeros((combined_df.shape[1], combined_df.shape[1])), 
                                      index=combined_df.columns, columns=combined_df.columns)
    p_value_matrix = pd.DataFrame(np.ones((combined_df.shape[1], combined_df.shape[1])), 
                                  index=combined_df.columns, columns=combined_df.columns)

    tasks = []
    with ProcessPoolExecutor() as executor:
        for i, col1 in enumerate(combined_df.columns):
            for j, col2 in enumerate(combined_df.columns):
                if i < j:
                    tasks.append(executor.submit(calculate_correlation_pair, col1, col2, combined_df))

        p_values = []
        index_pairs = []
        for future in as_completed(tasks):
            (col1, col2), corr, p_value = future.result()
            correlation_matrix.at[col1, col2] = corr
            correlation_matrix.at[col2, col1] = corr
            p_value_matrix.at[col1, col2] = p_value
            p_value_matrix.at[col2, col1] = p_value
            p_values.append(p_value)
            index_pairs.append((col1, col2))

    # Apply the chosen correction method
    reject, p_values_corrected, _, _ = multipletests(p_values, method=correction_method)

    for (col1, col2), p_value_corr in zip(index_pairs, p_values_corrected):
        p_value_matrix.at[col1, col2] = p_value_corr
        p_value_matrix.at[col2, col1] = p_value_corr

    correlation_matrix.columns = [col.replace(kme_prefix, '') for col in correlation_matrix.columns]
    correlation_matrix.index = [idx.replace(kme_prefix, '') for idx in correlation_matrix.index]

    p_value_matrix.columns = [col.replace(kme_prefix, '') for col in p_value_matrix.columns]
    p_value_matrix.index = [idx.replace(kme_prefix, '') for idx in p_value_matrix.index]

    return correlation_matrix, p_value_matrix, combined_df


def compare_plants(adatas: List[AnnData], correction_method='fdr_bh') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compares plants using aggregated eigengenes and returns the correlation matrix and p-value matrix.

    Parameters:
    - adatas (List[AnnData]): List of AnnData objects.
    - correction_method (str): Method to use for multiple testing correction.

    Returns:
    - Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: Correlation matrix, p-value matrix and module membership DF.
    """

    correlation_matrix, p_value_matrix, combined_df = calculate_correlation_matrix(
        adatas, correction_method=correction_method)

    return correlation_matrix, p_value_matrix, combined_df

# END CORRELATION OF EIGENGENES

# START HELPER FUNCTIONS FOR VISUALIZATION OF HOG LEVELS PER SPECIES


def extract_hog_levels_per_species(adatas: List[AnnData]) -> Dict[str, Dict[str, set]]:
    """
    Extracts and combines HOG levels from multiple AnnData objects per species.

    Parameters:
    - adatas (List[AnnData]): List of AnnData objects containing HOG information.

    Returns:
    - Dict[str, Dict[str, set]]: A dictionary where keys are species names and values are dictionaries 
            of HOG levels and their counts.
    """

    species_levels = {}
    for adata in adatas:
        species = adata.uns['species']
        if species not in species_levels:
            species_levels[species] = {}
        for col in adata.var.columns:
            if col.startswith('N') and col.endswith('_ortho'):
                level_name = col.replace('_ortho', '')
                if level_name not in species_levels[species]:
                    species_levels[species][level_name] = set()
                species_levels[species][level_name].update(
                    adata.var[col].dropna().unique())

    return species_levels


def get_common_hogs_per_level(species_hog_levels: Dict[str, Dict[str, set]]) -> Dict[str, set]:
    """
    Identifies common HOGs per level from the species HOG levels dictionary.

    Parameters:
    - species_hog_levels (Dict[str, Dict[str, set]]): A dictionary where keys are species names and values are dictionaries 
                                of HOG levels and their counts.

    Returns:
    - Dict[str, set]: A dictionary where keys are level names and values are sets of common HOGs.
    """

    common_hogs_per_level = {}
    levels = list(next(iter(species_hog_levels.values())).keys())

    for level in levels:
        common_hogs = set.intersection(
            *(species_hog_levels[species][level] for species in species_hog_levels))
        common_hogs_per_level[level] = common_hogs

    return common_hogs_per_level

# END HELPER FUNCTIONS FOR VISUALIZATION OF HOG LEVELS PER SPECIES

# START HOG EIGENGENES LINE PLOT BETWEEN TWO PLANTS AND MODULES


def aggregate_eigengenes_filter(adata: AnnData, 
                                module: str = 'module_colors', 
                                term: str = '', 
                                eigengenes_key: str = 'eigengenes', 
                                ortho_id_col: str = 'ortho_id', 
                                kme_prefix: str = 'kME') -> pd.DataFrame:
    """
    Aggregates eigengenes for the specified module/term.

    Parameters:
    - adata (AnnData): The AnnData object.
    - module (str): Column name in 'adata.var' to group by. 
    - term (str): Module term to filter by. 
    - eigengenes_key (str): Key for accessing eigengenes in adata.varm. 
    - ortho_id_col (str): Column name for ortholog IDs. 
    - kme_prefix (str): Prefix for kME columns. 

    Returns:
    - pd.DataFrame: Aggregated eigengenes DataFrame.
    """
    # Match eigengenes with ortho_ids and module terms
    eigengenes_df = pd.DataFrame(adata.varm[eigengenes_key], index=adata.var_names)
    eigengenes_df.columns = [col.replace(kme_prefix, '') for col in eigengenes_df.columns]
    eigengenes_df[ortho_id_col] = adata.var[ortho_id_col]
    eigengenes_df['module'] = adata.var[module]

    # Filter by the specified module term
    subset = eigengenes_df[eigengenes_df['module'] == term]
    if subset.empty:
        return pd.DataFrame()

    # Remove the 'module' column and aggregate by ortho_id
    subset = subset.drop(columns=['module'])
    aggregated_df = subset.groupby(ortho_id_col).mean()

    return aggregated_df


def get_common_hogs(aggregated_df1: pd.DataFrame, aggregated_df2: pd.DataFrame) -> pd.Index:
    """
    Returns the common HOGs between two DataFrames.

    Parameters:
    - aggregated_df1 (pd.DataFrame): First DataFrame of aggregated eigengenes.
    - aggregated_df2 (pd.DataFrame): Second DataFrame of aggregated eigengenes.

    Returns:
    - pd.Index: Common HOGs.
    """
    common_hogs = aggregated_df1.index.intersection(aggregated_df2.index)
    return common_hogs

# END HOG EIGENGENES LINE PLOT BETWEEN TWO PLANTS AND MODULES
