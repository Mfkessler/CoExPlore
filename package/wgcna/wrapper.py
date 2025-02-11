import os
from .ortho import transcript_ortho_browser, summarize_orthogroups, prepare_dataframe, calculate_jaccard_matrix
from .utils import (load_subset_from_hdf5, 
                    calculate_eigengenes, 
                    add_goea_to_anndata, 
                    get_goea_df, 
                    export_co_expression_network_to_cytoscape, 
                    identify_network_clusters_from_json, 
                    calculate_all_tom_metrics)
from .plotting import (PlotConfig, 
                       plot_co_expression_network, 
                       plot_module_trait_relationships_heatmap, 
                       plot_eigengenes, 
                       plot_filtered_go_terms, 
                       plot_cyto_network, 
                       plot_overlap, 
                       plot_eigengene_expression_bokeh)
from anndata import AnnData
from jinja2 import Environment, FileSystemLoader
from typing import Union, List, Tuple, Dict, Callable
import pandas as pd


def analyze_co_expression_network(adata: Union[AnnData, List[AnnData]], config: PlotConfig, custom_filename: str = None, transcripts: str = None, 
                                  query: Union[str, int, List[str], List[int], None] = "GO:0009908", 
                                  topic: str = "Flower development", threshold: float = 0.2,
                                  tom_prefix: str = "/vol/share/ranomics_app/data/tom", tom_path: str = None, goea_threshold: int = 10,
                                  node_threshold_percent: float = 0.02, plot_module_trait_relationships: bool = False,
                                  plot_bar_plots: bool = True, plot_go_enrichment: bool = True,
                                  out: str = "html", obo_path: str = "../go-basic.obo", 
                                  template_path: str = "../flask/app_dev/templates", highlight: List[str] = None,
                                  highlight_color: str = "magenta", tool: str = "cytoscape",
                                  use_colors: bool = False, use_shapes: bool = False,
                                  node_size: int = 10, use_symmetry: bool = False, progress_callback: Callable[[str], None] = None,
                                  trait: str = "tissue", filter_edges: bool = True, include_neighbours: bool = False,
                                  max_neighbors: int = 10) -> dict:
    """
    Analyze co-expression network for a given topic:
    - Plot the co-expression network and identify clusters.
    - Calculate eigengenes for clusters in the topic.
    - Plot module-trait relationships heatmap.
    - Plot bar plots of cluster eigengenes.
    - Perform GO enrichment analysis for clusters.
    - Generate an HTML report with the results.

    Parameters:
    - adata (Union[AnnData, List[AnnData]]): AnnData object(s) to analyze.
    - config (PlotConfig): Configuration object for plot control.
    - custom_filename (str): Custom filename for the output files.
    - transcripts (str): List of transcripts to analyze.
    - query (Union[str, int, List[str], List[int], None]): Query to search for transcripts.
    - topic (str): Headline for the analysis.
    - threshold (float): Threshold for the co-expression network (TOM value).
    - tom_prefix (str): Prefix for the TOM path.
    - goea_threshold (int): Threshold for the GO enrichment analysis.
    - node_threshold_percent (float): Percentage of nodes required in a cluster to be assigned.
    - plot_module_trait_relationships (bool): Plot module-trait relationships heatmap.
    - plot_bar_plots (bool): Plot bar plots of cluster eigengenes.
    - plot_go_enrichment (bool): Perform GO enrichment analysis for clusters.
    - out (str): Output type, can be "html" for HTML output.
    - obo_path (str): Path to the Gene Ontology OBO file.
    - template_path (str): Path to the Jinja2 template.
    - highlight (List[str]): List of transcripts to highlight in the co-expression network.
    - highlight_color (str): Color to use for highlighting nodes.
    - tool (str): Tool to use for plotting the co-expression network.
    - use_colors (bool): Use module colors for nodes of the co-expression network.
    - use_shapes (bool): Use shapes for modules in the co-expression network.
    - node_size (int): Size of nodes in the co-expression network.
    - use_symmetry (bool): Use one half of the TOM matrix to generate the network.
    - progress_callback (Callable[[str], None]): Callback function to report progress.
    - trait (str): Column name in adata.obs containing the trait information.
    - filter_edges (bool): Filter edges to increase performance.
    - include_neighbours (bool): If True, for each transcript in rows, 
                                check for neighbours with value >= threshold and include them.
    - max_neighbors (int): Maximum number of neighbours to include.

    Returns:
    - dict: Eigengenes for clusters if out is not "html".
    - str: Path to the generated HTML if out is "html".
    """

    if tool not in ["cytoscape", "plotly"]:
        raise ValueError(f"Invalid tool: {tool}, use 'cytoscape' or 'plotly'")
    
    # Check if adata is a list of AnnData objects
    if tool == "cytoscape" and isinstance(adata, list):
        if custom_filename is None:	
            custom_filename = f"{topic}"

    elif tool == "plotly" and isinstance(adata, list):
        raise ValueError("Plotly is not supported for multiple AnnData objects")
    else:
        if custom_filename is None:	
            custom_filename = f"{adata.uns['name']}_{topic}"

    print(f"Analyzing co-expression network for {topic} using {tool}")
    if progress_callback:
        progress_callback(f"Analyzing co-expression network for {topic} using {tool}")

    if isinstance(adata, list):
        tom, adata = get_tom_data(tom_path, adata, transcripts=transcripts, query=query, threshold=threshold, tom_prefix=tom_prefix, 
                                  use_symmetry=use_symmetry, progress_callback=progress_callback, include_neighbours=include_neighbours,
                                  max_neighbors=max_neighbors)
    else:
        tom = get_tom_data(tom_path, adata, transcripts=transcripts, query=query, threshold=threshold, tom_prefix=tom_prefix, 
                           use_symmetry=use_symmetry, progress_callback=progress_callback, include_neighbours=include_neighbours,
                           max_neighbors=max_neighbors)

    title_suffix = topic
    print(f"Plotting co-expression network for {topic}")
    if tool == "cytoscape":
        if progress_callback:
            progress_callback(f"Preparing data for {tool} visualization")
        cyto_data = export_co_expression_network_to_cytoscape(tom, adata, threshold=threshold)

        if progress_callback:
            progress_callback(f"Identifying clusters")
        cluster_map = identify_network_clusters_from_json(cyto_data, topic, node_threshold_percent=node_threshold_percent)

        ortho_df = prepare_dataframe(cyto_data['nodes'], cluster_map)
        ortho_table = summarize_orthogroups(ortho_df).reset_index()
        ortho_table_html = ortho_table.to_html(classes='dynamic-table display dataTable no-border', index=False, border=0, header=True, table_id="orthoResults")

        ortho_jaccaard_matrix = calculate_jaccard_matrix(ortho_df)
        plot_overlap(ortho_jaccaard_matrix, config, title="Jaccard Similarity Between Clusters Based on Orthogroups", 
                     x_label="Cluster", y_label="Cluster", custom_filename=f"{custom_filename}_overlap", height=None, width=None)

        if progress_callback:
            progress_callback(f"Plotting {tool} network")
        plot_cyto_network(config, custom_filename=custom_filename, network_data=cyto_data, use_colors=use_colors, use_shapes=use_shapes,
                          cluster_info=cluster_map, searchpath=template_path, node_size=node_size, highlight=highlight, 
                          highlight_color=highlight_color, filter_edges=filter_edges)
    else:
        if progress_callback:
            progress_callback(f"Plotting {tool} network")
        cluster_map = plot_co_expression_network(tom, adata, config, cluster_name=None, threshold=threshold, 
                                    use_shapes=use_shapes, use_colors=use_colors, title_suffix=title_suffix, 
                                    custom_filename=custom_filename, show_symbol_legend=True,
                                    show_edge_legend=True, node_size=node_size, identify_clusters=topic,
                                    node_threshold_percent=node_threshold_percent, height=None, width=None,
                                    highlight=highlight, highlight_color=highlight_color)
        
    # Calculate TOM metrics and generate Datatable
    if progress_callback:
        progress_callback(f"Calculating TOM metrics")
    tom_metrics = calculate_all_tom_metrics(tom, adata, cluster_map, tool, progress_callback).reset_index()
    tom_table_html = tom_metrics.to_html(classes='dynamic-table display dataTable no-border', index=False, border=0, header=True, table_id="tomResults")

    if progress_callback:
        progress_callback(f"Calculating eigengenes")
    eigengene_files, module_trait_files, table_html, combined_plot_path = process_eigengenes(adata, cluster_map, config, topic, 
                                                                          plot_module_trait_relationships, plot_bar_plots, 
                                                                          custom_filename=custom_filename, tool=tool, 
                                                                          progress_callback=progress_callback,
                                                                          trait=trait)

    if plot_go_enrichment:
        print(f"Performing GO enrichment analysis for clusters")
        add_goea_to_anndata(adata, column=topic, results_dir=f"{config.output_path}", uns_key=topic, top_percentage=goea_threshold, 
                            go_background_set=None, obo_path=obo_path)
        df = get_goea_df(adata, key=topic, column_name=f"{topic} Cluster")
        plot_filtered_go_terms(df, config, column=f"{topic} Cluster")

    if out == "html":
        # Prepare list of images (module_trait_file + eigengene_files)
        images = []
        if plot_module_trait_relationships:
            for png_file in module_trait_files:
                images.append(os.path.basename(png_file))
        
        if plot_bar_plots:
            for png_file in eigengene_files:
                images.append(os.path.basename(png_file))
        
        ortho_plot_path = f"{custom_filename}_overlap.html"
        network_plot_path = f"{custom_filename}.html"
        combined_plot_path = os.path.basename(combined_plot_path) if combined_plot_path else None
        
        if progress_callback:
            progress_callback(f"Generating HTML report")
        output_file_path = generate_co_expression_html(
            adata=adata,
            config=config,
            table_html=table_html,
            tom_table_html=tom_table_html,
            ortho_table_html=ortho_table_html if tool == "cytoscape" else None,
            images=images,
            combined_plot_path=combined_plot_path,
            topic=topic,
            ortho_plot_path=ortho_plot_path if tool == "cytoscape" else None,
            network_plot_path=network_plot_path,
            threshold=threshold,
            template_path=template_path,
            custom_filename=f"{custom_filename}_analysis.html"
        )
        
        return os.path.basename(output_file_path)
    else: 
        return "Not implemented"


def get_tom_data(tom_path: Union[str, List[str]], adata: Union[AnnData, List[AnnData]], transcripts: Dict[str, List[str]] = None, query: str = "GO:0009908", 
             threshold: float = 0.2, tom_prefix: str = "/vol/share/ranomics_app/data/tom", use_symmetry: bool = False, index_map: dict = None,
             columns_map: dict = None, progress_callback: Callable[[str], None] = None,
             include_neighbours: bool = False, max_neighbors: int = 10) -> Union[pd.DataFrame, List[pd.DataFrame]]:
    """
    This function loads the TOM matrix for one or multiple AnnData objects.
    
    Parameters:
    - tom_path (Union[str, List[str]]): Path to the TOM matrix file or a list of paths if multiple.
    - adata (Union[AnnData, List[AnnData]]): The AnnData object or a list of AnnData objects.
    - transcripts (Dict[str, List[str]]): List of transcripts to be used for subsetting.
    - query (str): The query for filtering transcripts. Defaults to "GO:0009908".
    - threshold (float): Threshold for filtering the TOM matrix. Defaults to 0.2.
    - tom_prefix (str): The prefix for the TOM file if tom_path is not provided.
    - use_symmetry (bool): Whether to use a symmetric TOM matrix.
    - index_map (dict): A dictionary mapping row labels to indices.
    - columns_map (dict): A dictionary mapping column labels to indices.
    - progress_callback (Callable[[str], None]): Callback function to report progress.
    - include_neighbours (bool): If True, for each transcript in rows,
                                check for neighbours with value >= threshold and include them.  
    - max_neighbors (int): Maximum number of neighbours to include.

    Returns:
    - tom (Union[np.ndarray, List[np.ndarray]]): The loaded TOM matrix or a list of TOM matrices.
    """

    valid_adatas = []

    if isinstance(adata, list):
        # Handle the case where adata is a list of AnnData objects
        toms = []
        for idx, ad in enumerate(adata):
            current_tom_path = tom_path[idx] if isinstance(tom_path, list) else f"{tom_prefix}/tom_matrix_{ad.uns['name']}.h5"

            if not os.path.exists(current_tom_path):
                print(f"File not found: {current_tom_path}")
                continue
            
            if not transcripts:
                current_transcripts = [t[1] for t in list(transcript_ortho_browser(query, ad).index)]

            if transcripts:
                current_species = ad.uns['species']
                if current_species in transcripts:
                    # Filter transcripts specific to the current species that are present in ad.var_names
                    current_transcripts = [t for t in transcripts[current_species] if t in ad.var_names]
                else:
                    print(f"No transcripts found for species {current_species} in dataset {ad.uns['name']}")
                    continue

            if not current_transcripts:
                print(f"No transcripts found for the query in dataset {ad.uns['name']}")
                continue
                
            print(f"Loading {len(current_transcripts)} transcripts from {current_tom_path}")
            if progress_callback:
                progress_callback(f"Loading {len(current_transcripts)} transcripts from TOM of dataset {ad.uns['name']}")
            tom = load_subset_from_hdf5(filename=current_tom_path, rows=current_transcripts, cols=current_transcripts, 
                                        threshold=threshold, use_symmetry=use_symmetry, index_map=index_map, columns_map=columns_map,
                                        include_neighbours=include_neighbours, max_neighbors=max_neighbors)
            print(f"Transcripts after filtering: {tom.shape[0]}")
            toms.append(tom)
            valid_adatas.append(ad)

        for idx, t in enumerate(toms):
            toms[idx] = t.astype("float64")

        return toms, valid_adatas
    else:
        if not tom_path:
            tom_path = f"{tom_prefix}/tom_matrix_{adata.uns['name']}.h5"

        if not os.path.exists(tom_path):
            print(f"File not found: {tom_path}")
            return None

        if transcripts:
            transcripts = [t for t in transcripts[adata.uns['species']] if t in adata.var_names]
        else:
            transcripts = [t[1] for t in list(transcript_ortho_browser(query, adata).index)]

        if not transcripts:
            print(f"No transcripts found for the query in dataset {adata.uns['name']}")
            return None
        
        print(f"Loading {len(transcripts)} transcripts from {tom_path}")
        if progress_callback:
            progress_callback(f"Loading {len(transcripts)} transcripts from TOM of dataset {adata.uns['name']}")
        tom = load_subset_from_hdf5(filename=tom_path, rows=transcripts, cols=transcripts, threshold=threshold, use_symmetry=use_symmetry,
                                    include_neighbours=include_neighbours, max_neighbors=max_neighbors)
        print(f"Transcripts after filtering: {tom.shape[0]}")

        tom = tom.astype("float64")

        return tom
    

def process_eigengenes(adata: Union[AnnData, List[AnnData]], cluster_map: dict, config: dict, topic: str, 
                       plot_module_trait_relationships: bool = False, plot_bar_plots: bool = False, 
                       custom_filename: str = None, tool: str = "cytoscape", progress_callback: Callable[[str], None] = None,
                       trait: str = "tissue") -> Tuple[List[str], List[str], str]:
    """
    This function calculates eigengenes for one or multiple AnnData objects and generates plots if requested.
    
    Parameters:
    - adata (Union[AnnData, List[AnnData]]): The AnnData object or a list of AnnData objects.
    - cluster_map (dict): A dictionary mapping transcripts to clusters.
    - config (dict): Configuration parameters for plotting.
    - topic (str): The column name for the topic used in eigengene calculation.
    - plot_module_trait_relationships (bool): Whether to plot module-trait relationships.
    - plot_bar_plots (bool): Whether to generate bar plots for eigengenes.
    - custom_filename (str, optional): A custom filename prefix for saving plots.
    - tool (str): The tool to use for plotting the co-expression network.
    - progress_callback (Callable[[str], None]): Callback function to report progress.
    - trait (str): Column name in adata.obs containing the trait information.

    Returns:
    - eigengene_files (List[str]): A list of paths to eigengene plot files.
    - module_trait_files (List[str]): A list of paths to module-trait relationship heatmaps.
    - table_html (str): The HTML string of the module information table.
    """

    eigengene_files = []
    module_trait_files = []
    all_module_info_dfs = []
    eigenene_data = []
    columns = ["Species", "Total Transcripts", "Top Module", "Top Module (%)", "Unique Modules"]
    combined_plot_path = None
    
    if isinstance(adata, list):
        # Handle the case where adata is a list of AnnData objects
        for ad in adata:
            if not trait in ad.obs.columns:
                if "Combined_Trait" not in ad.obs.columns:
                    raise KeyError(f"Trait column '{trait}' not found in adata.obs for dataset '{ad.uns['name']}'")
                else:
                    trait = "Combined_Trait"
            
            print(f"Calculating eigengenes for clusters in {topic} for dataset {ad.uns['name']}")
            if progress_callback:
                progress_callback(f"Calculating eigengenes for clusters of dataset {ad.uns['name']}")

            if tool == "cytoscape":
                # Filter cluster_map to only include clusters from the current species
                cluster_map_temp = {"_".join(k.split("_")[1:]): v for k, v in cluster_map.items() if k.startswith(ad.uns['species'])}

            # Check if the cluster_map is empty
            if not cluster_map_temp:
                print(f"No clusters found for dataset {ad.uns['name']}")
                print("Skipping eigengene calculation and plotting")
                continue

            eigengenes = calculate_eigengenes(ad, cluster_map_temp, column=topic)
            module_info_df = eigengenes['moduleInfo']
            cluster_eigengenes = eigengenes['eigengenes']

            if plot_module_trait_relationships:
                print(f"Plotting module-trait relationships heatmap for dataset {ad.uns['name']}")
                module_trait_file = plot_module_trait_relationships_heatmap(ad, cluster_eigengenes, module_info_df, config, 
                                                                            trait, topic, 
                                                                            custom_filename=f"{custom_filename}_{ad.uns['name']}")
                module_trait_files.append(module_trait_file)

            if plot_bar_plots:
                print(f"Plotting bar plots of cluster eigengenes for dataset {ad.uns['name']}")
                if trait not in ad.obs:
                    raise KeyError(f"Trait column '{trait}' not found in adata.obs for dataset '{ad.uns['name']}'")
                
                custom_stages_temp = ad.obs.get(trait, pd.Series()).unique().tolist() if "Combined_Trait" in ad.obs.columns else None

                if progress_callback:
                    progress_callback(f"Plotting bar plots of cluster eigengenes for dataset {ad.uns['name']}")
                eigengene_plot_files = plot_eigengenes(ad, cluster_eigengenes, config, column=topic, custom_filename=f"{custom_filename}_{ad.uns['name']}",
                                                       trait=trait, custom_stages=custom_stages_temp)
                eigengene_files.extend(eigengene_plot_files)

            module_info_df.columns = columns
            module_info_df = module_info_df.reset_index()  # Keep the index as a column (Cluster)
            all_module_info_dfs.append(module_info_df)

            # Convert sample index to column and rename
            cluster_eigengenes = map_tissue_info(cluster_eigengenes, ad)

            # Melt DataFrame to long format
            cluster_eigengenes = cluster_eigengenes.melt(id_vars=['Sample', 'Species', 'Tissue'], var_name='Cluster', value_name='Expression')

            # Append to the list
            eigenene_data.append(cluster_eigengenes)
    else:
        if not trait in adata.obs.columns:
            if "Combined_Trait" not in adata.obs.columns:
                raise KeyError(f"Trait column '{trait}' not found in adata.obs for dataset '{adata.uns['name']}'")
            else:
                trait = "Combined_Trait"
        # Handle the case where adata is a single AnnData object
        if tool == "cytoscape":
            # Filter cluster_map to only include clusters from the current species
            cluster_map = {k.split("_")[1]: v for k, v in cluster_map.items() if k.startswith(adata.uns['species'])}

        print(f"Calculating eigengenes for clusters in {topic}")
        if progress_callback:
            progress_callback(f"Calculating eigengenes for clusters of dataset {adata.uns['name']}")
        eigengenes = calculate_eigengenes(adata, cluster_map, column=topic)
        module_info_df = eigengenes['moduleInfo']
        cluster_eigengenes = eigengenes['eigengenes']

        if plot_module_trait_relationships:
            print(f"Plotting module-trait relationships heatmap")
            module_trait_file = plot_module_trait_relationships_heatmap(adata, cluster_eigengenes, module_info_df, config, 
                                                                        trait, topic, 
                                                                        custom_filename=custom_filename)
            module_trait_files.append(module_trait_file)

        if plot_bar_plots:
            print(f"Plotting bar plots of cluster eigengenes")

            if trait not in adata.obs:
                raise KeyError(f"Trait column '{trait}' not found in adata.obs for dataset '{adata.uns['name']}'")
            
            custom_stages_temp = adata.obs.get(trait, pd.Series()).unique().tolist() if "Combined_Trait" in adata.obs.columns else None

            if progress_callback:
                progress_callback(f"Plotting bar plots of cluster eigengenes")
            eigengene_files = plot_eigengenes(adata, cluster_eigengenes, config, column=topic, custom_filename=custom_filename,
                                              trait=trait, custom_stages=custom_stages_temp)

        module_info_df.columns = columns
        module_info_df = module_info_df.reset_index()  # Keep the index as a column (Cluster)
        all_module_info_dfs.append(module_info_df)

        cluster_eigengenes = map_tissue_info(cluster_eigengenes, adata)

        # Melt DataFrame to long format
        cluster_eigengenes = cluster_eigengenes.melt(id_vars=['Sample', 'Species', 'Tissue'], var_name='Cluster', value_name='Expression')

        # Append to the list
        eigenene_data.append(cluster_eigengenes)

    if progress_callback:
        progress_callback(f"Combining eigengene data for plotting")
    combined_eigengene_data = pd.concat(eigenene_data, ignore_index=True)
    if progress_callback:
        progress_callback(f"Plotting combined eigengene expression")
    combined_plot_path = plot_eigengene_expression_bokeh(combined_eigengene_data, config,
                                                    custom_filename=f"{custom_filename}_combined_bokeh")
    
    # Combine all module info DataFrames if needed
    combined_module_info_df = pd.concat(all_module_info_dfs, ignore_index=False)
    
    # Convert DataFrame to HTML table
    combined_module_info_df = combined_module_info_df.reset_index(drop=True)
    table_html = combined_module_info_df.to_html(classes='dynamic-table display dataTable no-border', index=False, border=0, header=True, table_id="infoResults")
    
    return eigengene_files, module_trait_files, table_html, combined_plot_path


def map_tissue_info(cluster_eigengenes: pd.DataFrame, ad: AnnData) -> pd.DataFrame:
    """
    Map tissue information from ad.obs to the 'cluster_eigengenes' DataFrame. 
    If 'tissue' column exists in ad.obs, use it. Otherwise, if 'Combined_Trait' 
    exists, use that. If none of these columns exist, raise an error.
    
    Parameters:
    - cluster_eigengenes (pd.DataFrame): DataFrame containing eigengenes with 'Sample' column.
    - ad (AnnData): The AnnData object containing metadata in ad.obs.
    
    Returns:
    - pd.DataFrame: The updated 'cluster_eigengenes' DataFrame with 'Tissue' column added.
    """

    # Rename column if needed
    cluster_eigengenes = cluster_eigengenes.reset_index()
    cluster_eigengenes = cluster_eigengenes.rename(columns={'sample': 'Sample'})

    # Check which column to use for Tissue information
    if 'tissue' in ad.obs.columns:
        cluster_eigengenes['Tissue'] = cluster_eigengenes['Sample'].map(ad.obs['tissue'])
    elif 'Combined_Trait' in ad.obs.columns:
        cluster_eigengenes['Tissue'] = cluster_eigengenes['Sample'].map(ad.obs['Combined_Trait'])
    else:
        raise ValueError("Neither 'tissue' nor 'Combined_Trait' columns found in ad.obs!")

    # Add Species from AnnData metadata
    cluster_eigengenes['Species'] = ad.uns['species']

    return cluster_eigengenes


def generate_co_expression_html(adata: Union[AnnData, List[AnnData]], config: PlotConfig, table_html: str, tom_table_html: str,
                                ortho_table_html: str, images: List[str], combined_plot_path: str, topic: str, ortho_plot_path: str, network_plot_path: str, 
                                threshold: int, template_path: str, custom_filename: str = None) -> str:
    """
    Generate an HTML report for the co-expression analysis.

    Parameters:
    - adata (Union[AnnData, List[AnnData]]): AnnData object to analyze.
    - config (PlotConfig): Configuration object for plot control.
    - table_html (str): HTML table with module information.
    - images (List[str]): List of image file names.
    - combined_plot_path (str): Path to the combined eigengene plot.
    - topic (str): Headline for the analysis.
    - ortho_plot_path (str): Path to the HTML file of the orthogroup overlap plot
    - network_plot_path (str): Path to the HTML file of the co-expression network plot.
    - threshold (float): Threshold for the co-expression network (TOM value).
    - template_path (str): Path to the Jinja2 template.
    - custom_filename (str): Custom filename for the output file.

    Returns:
    - str: Path to the generated HTML
    """

    if custom_filename is None:
        custom_filename = f"{topic}_{threshold}_analysis.html" if isinstance(adata, list) else f"{adata.uns['name']}_{topic}_{threshold}_analysis.html"

    env = Environment(loader=FileSystemLoader(searchpath=template_path))
    template = env.get_template("co_expr_template.html")
    
    rendered_html = template.render(
        topic=topic,
        table_html=table_html,
        tom_table_html=tom_table_html,
        ortho_table_html=ortho_table_html,
        images=images,
        combined_plot_path=combined_plot_path,
        network_plot_path=network_plot_path,
        ortho_plot_path=ortho_plot_path
    )
    
    output_file_path = os.path.join(config.output_path, custom_filename)
    
    with open(output_file_path, "w") as f:
        f.write(rendered_html)
    
    return output_file_path