import pandas as pd
import requests
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
    If delimiter is None, do NOT split value_col.
    """
    df = adata.var
    background_size = df.shape[0]

    df_exp = df[[group_col, value_col]].dropna()
    df_exp[value_col] = df_exp[value_col].astype(str)

    # Only split if delimiter is not None!
    if delimiter is not None:
        df_exp = df_exp.assign(**{value_col: df_exp[value_col].str.split(delimiter)})
        df_exp = df_exp.explode(value_col)
        df_exp[value_col] = df_exp[value_col].str.strip()
    else:
        # No splitting at all, value_col stays as-is
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
            proportion_in_group = a / grp_size if grp_size > 0 else 0
            results.append({
                group_col: grp,
                "term": term,
                "count_in_group": a,
                "count_in_total": count_in_all,
                "group_size": grp_size,
                "background_size": background_size,
                "proportion_in_group": proportion_in_group,
                "pvalue": pval
            })
    res_df = pd.DataFrame(results)
    if not res_df.empty:
        res_df['p_adj_bonf'] = res_df['pvalue'] * res_df.shape[0]
        res_df = res_df.sort_values(['pvalue', 'count_in_group'], ascending=[True, False])
        res_df = res_df[res_df['pvalue'] <= pval_thresh]

    return res_df.reset_index(drop=True)


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
    matrix_type: str = 'proportion',
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
            - 'delimiter': delimiter for multiple entries (e.g. ',' or '|')
            - 'mapping_func': function to map term->metadata
            - 'mapping_columns': columns for mapped values
            - 'prefix': prefix for feature names (else value_col is used)
    group_col : str
        Column to group by (default: 'module_labels')
    matrix_type : str
        Type of matrix to build: 'proportion', 'membership', or 'score'.
        - 'proportion': Proportion of terms in group
        - 'membership': Binary membership matrix (1 if term in group, else 0)
        - 'score': Score matrix with -log10(p_adj_bonf) for enrichment
    sort_cols_base : list of str
        Default sorting columns for enrichment results (used if not given per config)
    top_n : int
        How many top entries per group/category?

    Returns:
    AnnData : Feature matrix across all annotation types, modules as rows, features as columns.
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
            if 'go_cat_code' in enrich.columns:
                groupby_cols.append('go_cat_code')
            top_terms = enrich.groupby(groupby_cols).head(top_n).reset_index(drop=True)
            annots_dict[name] = top_terms
        all_annots.append(annots_dict)
        all_prefixes.append(prefix)

    # 2. Build feature matrices for each annotation type, depending on matrix_type
    def collect_unique_terms(annot_dict, value_col='term'):
        all_terms = set()
        for df in annot_dict.values():
            if value_col in df.columns:
                all_terms.update(df[value_col].dropna().unique())
        return sorted(all_terms)

    def build_score_matrix_adj(annot_dict, group_col, value_col, score_col='p_adj_bonf', all_values=None, prefix="", pval_thresh=0.05):
        if all_values is None:
            all_values = collect_unique_terms(annot_dict, value_col)
        module_names = []
        records = []
        scores = []
        for species, df in annot_dict.items():
            for _, row in df.iterrows():
                module = f"{species}_{row[group_col]}"
                module_names.append(module)
                records.append((module, row[value_col]))
                scores.append(row.get(score_col, np.nan))
        module_names = sorted(set(module_names))
        matrix = pd.DataFrame(0.0, index=module_names, columns=[f"{prefix}{col}" for col in all_values])
        for (module, term), score in zip(records, scores):
            col = f"{prefix}{term}"
            if col in matrix.columns:
                if pd.notnull(score) and score <= pval_thresh and score > 0:
                    val = -np.log10(score)
                    matrix.loc[module, col] = val
                else:
                    matrix.loc[module, col] = 0.0
        return matrix

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

    def build_proportion_matrix(annot_dict, group_col, value_col, score_col='proportion_in_group', all_values=None, prefix=""):
        if all_values is None:
            all_values = collect_unique_terms(annot_dict, value_col)
        module_names = []
        records = []
        scores = []
        for species, df in annot_dict.items():
            for _, row in df.iterrows():
                module = f"{species}_{row[group_col]}"
                module_names.append(module)
                records.append((module, row[value_col]))
                scores.append(row.get(score_col, 0))
        module_names = sorted(set(module_names))
        matrix = pd.DataFrame(0.0, index=module_names, columns=[f"{prefix}{col}" for col in all_values])
        for (module, term), score in zip(records, scores):
            col = f"{prefix}{term}"
            if col in matrix.columns:
                matrix.loc[module, col] = score  # Value between 0 and 1
        return matrix

    def build_combined_matrix(
        annot_dict, group_col, value_col,
        prop_col='proportion_in_group', score_col='p_adj_bonf',
        all_values=None, prefix="", pval_thresh=0.05
    ):
        if all_values is None:
            all_values = collect_unique_terms(annot_dict, value_col)
        module_names = []
        records = []
        props = []
        scores = []
        for species, df in annot_dict.items():
            for _, row in df.iterrows():
                module = f"{species}_{row[group_col]}"
                module_names.append(module)
                records.append((module, row[value_col]))
                props.append(row.get(prop_col, 0))
                s = row.get(score_col, np.nan)
                # -log10(p) nur wenn signifikant und p>0, sonst 0
                score = -np.log10(s) if pd.notnull(s) and s <= pval_thresh and s > 0 else 0.0
                scores.append(score)
        module_names = sorted(set(module_names))
        matrix = pd.DataFrame(0.0, index=module_names, columns=[f"{prefix}{col}" for col in all_values])
        for (module, term), prop, score in zip(records, props, scores):
            col = f"{prefix}{term}"
            if col in matrix.columns:
                matrix.loc[module, col] = prop * score
        return matrix

    # Build feature matrices according to matrix_type
    feature_matrices = []
    for annots_dict, config, prefix in zip(all_annots, config_list, all_prefixes):
        value_col = config['value_col']
        if matrix_type == "proportion":
            M = build_proportion_matrix(annots_dict, group_col=group_col, value_col='term', prefix=prefix)
        elif matrix_type == "membership":
            M = build_membership_matrix(annots_dict, group_col=group_col, value_col='term', prefix=prefix)
        elif matrix_type == "score":
            M = build_score_matrix_adj(annots_dict, group_col=group_col, value_col='term', prefix=prefix)
        elif matrix_type == "combined":
            M = build_combined_matrix(annots_dict, group_col=group_col, value_col='term', prefix=prefix)
        else:
            raise ValueError(f"Unknown matrix_type: {matrix_type!r}")
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
