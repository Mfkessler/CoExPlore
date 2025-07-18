import pandas as pd
import requests
import re
import numpy as np
from sklearn.decomposition import PCA
from umap import UMAP
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from scipy.stats import fisher_exact
import anndata as ad
from typing import List, Optional, Any, Tuple, Dict


def split_column(series: pd.Series, delimiter: Optional[str]) -> pd.Series:
    """
    Split a pandas Series by the specified delimiter.

    Parameters:
    - series (pd.Series): Input pandas Series to split.
    - delimiter (Optional[str]): Delimiter string to use for splitting. If None, returns original series.

    Returns:
    - pd.Series: Series with lists as values if delimiter was applied, otherwise original series.
    """

    if delimiter is None:
        return series
    if series.apply(lambda x: isinstance(x, list)).any():
        return series
    return series.astype(str).str.split(delimiter)


def extract_top_n_by_group(
    df: pd.DataFrame,
    group_col: str,
    value_col: str,
    n: int = 5,
    delimiter: Optional[str] = ','
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extract the top N most frequent items in value_col per group_col.

    Parameters:
    - df (pd.DataFrame): Input DataFrame containing the data.
    - group_col (str): Column to group by.
    - value_col (str): Column with values to count frequency.
    - n (int): Number of top items per group to extract.
    - delimiter (Optional[str]): Delimiter for splitting values. If None, assumes no splitting needed.

    Returns:
    - Tuple[pd.DataFrame, pd.DataFrame]: Tuple containing the top-N DataFrame and a pivoted summary table.
    """

    prepared = df[[group_col, value_col]].dropna(subset=[value_col]).copy()
    prepared[value_col] = split_column(prepared[value_col], delimiter)
    prepared = prepared.explode(value_col)
    value_counts = (
        prepared.groupby([group_col, value_col])
        .size()
        .reset_index(name='count')
    )
    topn_per_group = (
        value_counts.sort_values([group_col, 'count'], ascending=[True, False])
        .groupby(group_col)
        .head(n)
        .reset_index(drop=True)
    )
    pivot = topn_per_group.pivot_table(
        index=group_col,
        columns=topn_per_group.groupby(group_col).cumcount() + 1,
        values=[value_col, 'count'],
        aggfunc='first'
    )
    pivot.columns = [f'Top{num}_{val}' for val, num in pivot.columns]
    pivot.reset_index(inplace=True)
    return topn_per_group, pivot


def analyze_uniqueness_and_overlap(
    topn_per_group: pd.DataFrame,
    group_col: str,
    value_col: str
) -> pd.DataFrame:
    """
    Analyze the uniqueness and overlap of top-N entries in each group.

    Parameters:
    - topn_per_group (pd.DataFrame): DataFrame with top-N values per group.
    - group_col (str): Column containing group identifiers.
    - value_col (str): Column containing values to check for uniqueness/overlap.

    Returns:
    - pd.DataFrame: DataFrame with columns for group, value, uniqueness, and overlap count.
    """

    group_to_topn = (
        topn_per_group.groupby(group_col)[value_col]
        .apply(set)
        .to_dict()
    )
    unique_per_group = {
        grp: {val for val in vals if not any(
            val in other for g, other in group_to_topn.items() if g != grp)}
        for grp, vals in group_to_topn.items()
    }
    overlap_count = {
        grp: {val: sum(val in other for g, other in group_to_topn.items() if g != grp)
              for val in vals}
        for grp, vals in group_to_topn.items()
    }
    overlap_df = pd.DataFrame([
        {
            group_col: grp,
            'value': val,
            'unique': val in unique_per_group[grp],
            'overlap_count': overlap_count[grp][val]
        }
        for grp in group_to_topn
        for val in group_to_topn[grp]
    ])
    return overlap_df


def fetch_go_names(
    go_ids: List[str],
    requests_module: Any = requests,
    batch_size: int = 100
) -> Dict[str, str]:
    """
    Fetch GO term names for a list of GO IDs using the QuickGO API.

    Parameters:
    - go_ids (List[str]): List of GO term IDs (e.g., ["GO:0008150"]).
    - requests_module (Any): Requests module or compatible object for HTTP requests.
    - batch_size (int): Maximum number of IDs per API request.

    Returns:
    - Dict[str, str]: Mapping from GO ID to GO term name.
    """

    go_name_dict = {}
    for i in range(0, len(go_ids), batch_size):
        batch_ids = go_ids[i:i+batch_size]
        joined_ids = ",".join(batch_ids)
        url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{joined_ids}"
        resp = requests_module.get(url)
        if resp.status_code == 200:
            for entry in resp.json()["results"]:
                go_id = entry["id"]
                go_name = entry["name"]
                go_name_dict[go_id] = go_name
        else:
            print("Error with QuickGO request:", resp.status_code, resp.text)
    return go_name_dict


def id_type(val: str) -> str:
    """
    Determine the type of biological identifier (GO, InterPro, Pfam, or other).

    Parameters:
    - val (str): Identifier string to test.

    Returns:
    - str: One of 'go', 'interpro', 'pfam', or 'other'.
    """

    if isinstance(val, str):
        if val.startswith("GO:"):
            return "go"
        elif re.match(r"^IPR\d+$", val):
            return "interpro"
        elif re.match(r"^PF\d+$", val):
            return "pfam"
    return "other"


def summarize_top_n_features(
    adata: Any,
    group_col: str = 'module_labels',
    value_col: str = 'go_terms',
    n: int = 5,
    delimiter: Optional[str] = ',',
    id_to_desc: Optional[dict] = None,
    requests_module: Any = requests,
) -> pd.DataFrame:
    """
    Extract and annotate the top-N most frequent features per group from an AnnData object.

    If value_col is 'go_terms', fetches feature names from QuickGO.
    If value_col is 'ipr_id' or 'pfam_id', uses the corresponding desc column in adata.var as mapping,
    unless a custom id_to_desc mapping is provided.

    Parameters:
    - adata (Any): AnnData object with .var attribute as DataFrame.
    - group_col (str): Column to group by.
    - value_col (str): Column containing feature or annotation IDs.
    - n (int): Number of top features per group.
    - delimiter (Optional[str]): Delimiter to split value_col entries.
    - id_to_desc (Optional[dict]): Optional custom mapping from ID to description.
    - requests_module (Any): Requests module or compatible object for HTTP requests.

    Returns:
    - pd.DataFrame: DataFrame summarizing top-N features with annotation and group size.
    """
    print("Species:", adata.uns.get('species', 'N/A'))
    df = adata.var

    topn_df, pivot = extract_top_n_by_group(
        df, group_col=group_col, value_col=value_col, n=n, delimiter=delimiter
    )

    overlap_df = analyze_uniqueness_and_overlap(
        topn_df, group_col=group_col, value_col=value_col
    )
    overlap_df_sorted = overlap_df.sort_values(
        [group_col, 'unique', 'overlap_count'],
        ascending=[True, False, True]
    )
    min_overlap_per_group = overlap_df_sorted.groupby(
        group_col).first().reset_index()

    # Feature name handling
    if value_col == "go_terms":
        id_types = min_overlap_per_group['value'].map(id_type)
        name_dict = {}
        values = min_overlap_per_group['value'].unique()
        if 'go' in id_types.values:
            go_ids = [val for val in values if id_type(val) == 'go']
            name_dict.update(fetch_go_names(go_ids, requests_module=requests_module))
        min_overlap_per_group['feature_name'] = min_overlap_per_group['value'].map(name_dict).astype('object')

    elif value_col in ["ipr_id", "pfam_id"]:
        # Use id_to_desc if given, else auto-extract from adata.var
        if id_to_desc is None:
            # Find desc column name
            desc_col = "ipr_desc" if value_col == "ipr_id" else "pfam_desc"
            id_to_desc = {}
            for idx, row in df.iterrows():
                ids = str(row[value_col]).split(",")
                descs = str(row[desc_col]).split("|")
                for id_, desc_ in zip(ids, descs):
                    id_ = id_.strip()
                    desc_ = desc_.strip()
                    if id_ and desc_:
                        id_to_desc[id_] = desc_
        # Map each value (may be multiple, comma-separated) to its description(s)
        def map_ids_to_descs(ids: str) -> str:
            if not isinstance(ids, str) or not ids.strip():
                return ""
            descs = [id_to_desc.get(i.strip(), "") for i in ids.split(",")]
            return "|".join([d for d in descs if d])
        min_overlap_per_group['feature_name'] = min_overlap_per_group['value'].apply(map_ids_to_descs)

    else:
        min_overlap_per_group['feature_name'] = None

    # Add group size
    group_size = df[group_col].value_counts().to_dict()
    min_overlap_per_group['group_size'] = min_overlap_per_group[group_col].map(group_size)

    return min_overlap_per_group


def annotate_modules_with_weight(
    go_df: pd.DataFrame,
    ipr_df: pd.DataFrame,
    pfam_df: pd.DataFrame,
    module_col: str = "module_labels",
    feature_col: str = "feature_name",
    unique_col: str = "unique"
) -> pd.DataFrame:
    """
    Create a summary annotation for each module based on unique GO, InterPro, and Pfam annotations,
    and compute a weight for annotation confidence.
    
    Parameters:
    - go_df, ipr_df, pfam_df: DataFrames with module-level results, as produced by summarize_top_n_features.
    - module_col (str): Column indicating the module label.
    - feature_col (str): Column with the feature/annotation name.
    - unique_col (str): Column indicating whether annotation is unique in the module.

    Returns:
    - pd.DataFrame: DataFrame with columns ['Module', 'Annotation', 'Weight'].
    """
    modules = sorted(set(go_df[module_col]) | set(ipr_df[module_col]) | set(pfam_df[module_col]))
    records = []
    
    for mod in modules:
        # Try to find unique annotation per source
        go_row = go_df[(go_df[module_col] == mod) & (go_df[unique_col] == True)]
        ipr_row = ipr_df[(ipr_df[module_col] == mod) & (ipr_df[unique_col] == True)]
        pfam_row = pfam_df[(pfam_df[module_col] == mod) & (pfam_df[unique_col] == True)]
        
        sources = {
            "GO": go_row[feature_col].iloc[0] if not go_row.empty else None,
            "IPR": ipr_row[feature_col].iloc[0] if not ipr_row.empty else None,
            "PFAM": pfam_row[feature_col].iloc[0] if not pfam_row.empty else None,
        }
        n_unique = sum(1 for v in sources.values() if v is not None)
        
        # Fallback: if nothing unique, use the most frequent/first from each source
        if n_unique == 0:
            go_row = go_df[go_df[module_col] == mod]
            ipr_row = ipr_df[ipr_df[module_col] == mod]
            pfam_row = pfam_df[pfam_df[module_col] == mod]
            sources = {
                "GO": go_row[feature_col].iloc[0] if not go_row.empty else None,
                "IPR": ipr_row[feature_col].iloc[0] if not ipr_row.empty else None,
                "PFAM": pfam_row[feature_col].iloc[0] if not pfam_row.empty else None,
            }
        
        # Score
        if n_unique == 3:
            weight = 1.0
        elif n_unique == 2:
            weight = 2/3
        elif n_unique == 1:
            weight = 1/3
        else:
            weight = 0.0
        
        # Join for annotation (with label)
        annotation = "; ".join([f"{k}: {v}" for k, v in sources.items() if v])
        if not annotation:
            annotation = "uncharacterized"
        
        records.append({"Module": mod, "Annotation": annotation, "Weight": weight})

    return pd.DataFrame.from_records(records)


def generic_enrichment(
    adata: ad.AnnData,
    group_col: str,
    value_col: str,
    delimiter: str = ',',
    min_count: int = 2,
    pval_thresh: float = 0.05
) -> pd.DataFrame:
    """
    Perform enrichment analysis for annotation terms (GO, IPR, PFAM, etc.) in each group of an AnnData object.
    
    Parameters:
    - adata (AnnData): AnnData object with .var DataFrame containing group_col and value_col.
    - group_col (str): Column name in .var for grouping (e.g., 'module_labels').
    - value_col (str): Column name in .var with term lists (comma- or delimiter-separated).
    - delimiter (str): Delimiter for splitting term strings.
    - min_count (int): Minimum number of occurrences in group to test.
    - pval_thresh (float): P-value cutoff for reporting enrichment.
    
    Returns:
    - pd.DataFrame: Enriched terms per group with counts and p-values.
    """
    df = adata.var
    background_size = df.shape[0]

    # Explode the value_col so each row is a single term
    df_exp = df[[group_col, value_col]].dropna()
    df_exp[value_col] = df_exp[value_col].astype(str)
    df_exp = df_exp.assign(**{value_col: df_exp[value_col].str.split(delimiter)})
    df_exp = df_exp.explode(value_col)
    df_exp[value_col] = df_exp[value_col].str.strip()
    df_exp = df_exp[df_exp[value_col] != '']

    all_terms = df_exp[value_col].tolist()
    counter_all = Counter(all_terms)

    results = []
    for grp, sub in df_exp.groupby(group_col):
        terms_in_grp = sub[value_col].tolist()
        grp_size = df[df[group_col] == grp].shape[0]
        counter_grp = Counter(terms_in_grp)
        for term, count_in_grp in counter_grp.items():
            if count_in_grp < min_count:
                continue
            count_in_all = counter_all[term]
            a = count_in_grp  # term in group
            b = count_in_all - a  # term in rest
            c = grp_size - a      # no term in group
            d = background_size - a - b - c  # no term in rest
            table = [[a, b], [c, d]]
            _, pval = fisher_exact(table, alternative='greater')
            results.append({
                group_col: grp,
                "term": term,
                "count_in_group": a,
                "count_in_total": count_in_all,
                "group_size": grp_size,
                "background_size": background_size,
                "pvalue": pval
            })
    res_df = pd.DataFrame(results)
    if not res_df.empty:
        res_df['p_adj_bonf'] = res_df['pvalue'] * res_df.shape[0]
        res_df = res_df.sort_values(['pvalue', 'count_in_group'], ascending=[True, False])
        res_df = res_df[res_df['pvalue'] <= pval_thresh]
    return res_df.reset_index(drop=True)


def top_enriched_term_per_group(
    enrich_df: pd.DataFrame,
    group_col: str = "module_labels",
    term_col: str = "term"
) -> pd.DataFrame:
    """
    Reduce enrichment result to the top (most significantly enriched) term per group.

    Parameters:
    - enrich_df (pd.DataFrame): Output from generic_enrichment().
    - group_col (str): The column indicating groups/modules.
    - term_col (str): The column with the term/ID.

    Returns:
    - pd.DataFrame: DataFrame with the top term for each group/module.
    """
    if enrich_df.empty:
        return enrich_df

    # Sort: lowest p-value first, then highest count_in_group
    enrich_df = enrich_df.sort_values(
        ["pvalue", "count_in_group"],
        ascending=[True, False]
    )
    # Keep only the top term per group
    top_df = enrich_df.groupby(group_col).first().reset_index()
    return top_df[[group_col, term_col, "pvalue", "count_in_group", "group_size"]]


def fetch_go_names_and_categories(
    go_ids: List[str],
    batch_size: int = 100
) -> Dict[str, Tuple[str, Optional[str], Optional[str]]]:
    """
    Fetch GO term names and categories for a list of GO IDs using the QuickGO API.

    Parameters:
    - go_ids (List[str]): List of GO term IDs (e.g., ["GO:0008150"]).
    - batch_size (int): Maximum number of IDs per API request.

    Returns:
    - Dict[str, Tuple[str, Optional[str], Optional[str]]]: 
        Mapping from GO ID to a tuple of (GO term name, category code, category name).
        Category code is one of "P" (biological_process), "F" (molecular_function), "C" (cellular_component).
    """
    go_dict: Dict[str, Tuple[str, Optional[str], Optional[str]]] = {}
    for i in range(0, len(go_ids), batch_size):
        batch = go_ids[i:i+batch_size]
        url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/" + ",".join(batch)
        resp = requests.get(url)
        if resp.status_code == 200:
            for entry in resp.json().get("results", []):
                go_id = entry.get("id")
                go_name = entry.get("name")
                cat_name = entry.get("aspect")
                cat_code = cat_name[:1].upper() if cat_name else None  # "P", "F", "C"
                if go_id and go_name and cat_name:
                    go_dict[go_id] = (go_name, cat_code, cat_name)
        else:
            print(f"Error fetching GO info: {resp.status_code} {resp.text}")
    return go_dict


def collect_module_enrichment(
    adatas: List[ad.AnnData],
    config_list: List[Dict[str, Any]],
    group_col: str = 'module_labels',
    sort_cols_base: Optional[List[str]] = None,
    top_n: int = 1
) -> ad.AnnData:
    """
    Automated pipeline for module enrichment and feature matrix construction.
    
    Parameters:
    adatas : list
        List of AnnData objects.
    config_list : list of dicts
        Each dict describes one annotation type. Required keys:
            - 'value_col': column in .var or .obs with annotation (e.g. 'go_terms')
            - 'delimiter': (optional) delimiter for multiple entries (e.g. ',' or '|')
            - 'mapping_func': (optional) function to map term->metadata
            - 'mapping_columns': (optional) columns for mapped values
            - 'prefix': (optional) prefix for feature names (else value_col is used)
    group_col : str
        Column to group by (default: 'module_labels')
    sort_cols_base : list of str
        Default sorting columns for enrichment results (used if not given per config)
    top_n : int
        How many top entries per group/category?
    
    Returns:
    AnnData : Feature matrix across all annotation types, modules as rows, binary features.
    """

    if sort_cols_base is None:
        sort_cols_base = [group_col, 'pvalue', 'count_in_group']

    all_annots = []
    all_prefixes = []

    # 1. Collect enrichment results for each annotation type
    for config in config_list:
        value_col = config['value_col']
        delimiter = config.get('delimiter')
        mapping_func = config.get('mapping_func')
        mapping_columns = config.get('mapping_columns')
        sort_cols = config.get('sort_cols', sort_cols_base)
        prefix = config.get('prefix', value_col.upper() + "_")
        
        annots_dict = {}
        for adata in adatas:
            name = adata.uns.get('name', 'unknown')
            enrich = generic_enrichment(
                adata, group_col=group_col, value_col=value_col, delimiter=delimiter
            )
            if mapping_func is not None and mapping_columns is not None:
                unique_terms = [t for t in enrich['term'].unique() if isinstance(t, str)]
                term_info = mapping_func(unique_terms)
                for i, col in enumerate(mapping_columns):
                    enrich[col] = enrich['term'].map(lambda x: term_info.get(x, (None,) * len(mapping_columns))[i])
            # Sort and take top N per group/category
            enrich = enrich.sort_values(sort_cols, ascending=[True]*len(sort_cols))
            groupby_cols = [group_col]
            if 'go_cat_code' in enrich.columns:  # For GO only, but generalizes
                groupby_cols.append('go_cat_code')
            top_terms = enrich.groupby(groupby_cols).head(top_n).reset_index(drop=True)
            annots_dict[name] = top_terms
        all_annots.append(annots_dict)
        all_prefixes.append(prefix)

    # 2. Build feature matrices (membership matrices) for each annotation type
    feature_matrices = []
    def collect_unique_terms(annot_dict, value_col='term'):
        all_terms = set()
        for df in annot_dict.values():
            if value_col in df.columns:
                all_terms.update(df[value_col].dropna().unique())
        return sorted(all_terms)
    
    def build_membership_matrix(annot_dict, group_col, value_col, all_values=None, prefix=""):
        if all_values is None:
            all_values = collect_unique_terms(annot_dict, value_col)
        module_names = []
        records = []
        for species, df in annot_dict.items():
            for _, row in df.iterrows():
                module = f"{species}_{row[group_col]}"
                module_names.append(module)
                records.append((module, row[value_col]))
        module_names = sorted(set(module_names))
        matrix = pd.DataFrame(0, index=module_names, columns=[f"{prefix}{col}" for col in all_values])
        for module, term in records:
            col = f"{prefix}{term}"
            if col in matrix.columns:
                matrix.loc[module, col] = 1
        return matrix

    for annots_dict, config, prefix in zip(all_annots, config_list, all_prefixes):
        value_col = config['value_col']
        M = build_membership_matrix(annots_dict, group_col=group_col, value_col='term', prefix=prefix)
        feature_matrices.append(M)

    # 3. Intersect index for all feature matrices (modules present in all)
    common_index = feature_matrices[0].index
    for M in feature_matrices[1:]:
        common_index = common_index.intersection(M.index)
    feature_matrices = [M.loc[common_index] for M in feature_matrices]
    
    # 4. Concatenate all feature matrices column-wise
    M_combined = pd.concat(feature_matrices, axis=1)

    # 5. Convert to AnnData
    adata_features = ad.AnnData(X=np.array(M_combined.values, dtype=np.float32))
    adata_features.obs_names = M_combined.index.astype(str)
    adata_features.var_names = M_combined.columns.astype(str)

    return adata_features


def module_clustering_analysis(
    adata: ad.AnnData,
    eps: float = 0.7,
    min_samples: int = 3,
    pca_components: int = 2,
    umap_components: int = 2,
    umap_random_state: int = 0,
    figsize: Tuple[int, int] = (7, 5)
) -> Dict[str, Any]:
    """
    Perform dimensionality reduction (PCA, UMAP), clustering (DBSCAN) and visualization for module-level features.

    Parameters:
    adata : ad.AnnData
        AnnData object with modules as observations (rows) and features (columns).
    eps : float
        DBSCAN epsilon parameter (neighborhood radius).
    min_samples : int
        DBSCAN min_samples parameter (minimum samples for core point).
    pca_components : int
        Number of components for PCA.
    umap_components : int
        Number of components for UMAP.
    umap_random_state : int
        Random state for UMAP reproducibility.
    figsize : tuple
        Figure size for the cluster plot.

    Returns:
    results : dict
        Dictionary with the following keys:
            - 'pca': PCA projection (np.ndarray)
            - 'umap': UMAP projection (np.ndarray)
            - 'clusters': DBSCAN cluster labels (np.ndarray)
            - 'df_cluster': DataFrame with module, species, and cluster
    """
    # PCA
    pca = PCA(n_components=pca_components)
    proj_pca = pca.fit_transform(adata.X)

    # UMAP
    umap = UMAP(n_components=umap_components, random_state=umap_random_state)
    proj_umap = umap.fit_transform(adata.X)

    # Module and species extraction
    module_names = adata.obs_names
    species_per_module = [str(name).split("_")[0] for name in module_names]
    unique_species = sorted(set(species_per_module))

    palette = sns.color_palette("tab10", len(unique_species))
    species_colors = [palette[unique_species.index(s)] for s in species_per_module]

    # DBSCAN clustering
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    clusters = dbscan.fit_predict(proj_umap)

    # Plotting
    plt.figure(figsize=figsize)
    plt.scatter(proj_umap[:, 0], proj_umap[:, 1], c=clusters, cmap="tab20", s=50)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title("Module-Clustering (DBSCAN on UMAP)")
    plt.colorbar(label="Cluster")
    plt.show()

    # Cluster DataFrame
    df_cluster = pd.DataFrame({
        "module": module_names,
        "species": species_per_module,
        "cluster": clusters
    })

    print("\nModules per cluster:")
    print(df_cluster['cluster'].value_counts().sort_index())

    print("\nSingletons (Cluster -1):")
    print(df_cluster[df_cluster['cluster'] == -1][['module', 'species']])

    print("\nConserved clusters (with modules from multiple species):")
    for c in sorted(df_cluster['cluster'].unique()):
        if c == -1: continue
        subset = df_cluster[df_cluster['cluster'] == c]
        arts = set(subset['species'])
        if len(arts) > 1:
            print(f"  Cluster {c}: {len(subset)} modules from species {sorted(arts)}")

    return {
        "pca": proj_pca,
        "umap": proj_umap,
        "clusters": clusters,
        "df_cluster": df_cluster
    }
