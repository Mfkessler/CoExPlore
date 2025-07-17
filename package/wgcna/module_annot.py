import pandas as pd
import requests
import re
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
