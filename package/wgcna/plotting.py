import os
import networkx as nx
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scanpy as sc
import scipy.stats as stats
import json
from typing import List, Dict, Tuple
from anndata import AnnData
from .utils import reduce_tom_matrix, load_subset_from_hdf5, identify_network_clusters, generate_stage_color_dict, get_data_trait_df
from .ortho import extract_hog_levels_per_species, get_common_hogs_per_level, get_common_hogs, aggregate_eigengenes_filter
from plotly.subplots import make_subplots
from plotly.graph_objs import Figure
from matplotlib_venn import venn3
from matplotlib.colors import to_rgba
from upsetplot import UpSet
from sklearn.decomposition import PCA
from plotly.colors import sample_colorscale
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.godag_plot import plot_gos
from typing import List, Dict, Union
from scipy.stats import zscore
from jinja2 import Environment, FileSystemLoader
from scipy.cluster.hierarchy import linkage, leaves_list
from IPython.display import display, IFrame
from bokeh.plotting import figure, save, show, output_file
from bokeh.models import (
    ColumnDataSource, FactorRange, HoverTool, CustomJS,
    MultiSelect, Select, Whisker
)
from bokeh.layouts import column
from bokeh.transform import dodge, factor_cmap
from bokeh.palettes import Category20, viridis


plt.rcParams['savefig.bbox'] = 'tight'


class PlotConfig:
    def __init__(self, save_plots=False, save_raw=False, output_path='./', file_format='png', show=True, dpi=300):
        self.save_plots = save_plots
        self.save_raw = save_raw
        self.output_path = output_path
        self.file_format = file_format
        self.show = show
        self.dpi = dpi


###########
"""
WGCNA OBJECT PLOTS
Plotting functions for WGCNA objects using large matrices like TOM.
"""
###########


def plot_static_TOM_network(tom_path: str, adata: AnnData, config: PlotConfig, transcripts: List[str], selected_module_colors: List[str] = None,
                            threshold: float = 0.5, custom_filename: str = None) -> None:
    """
    Plots the co-expression network for large datasets using NetworkX and Matplotlib.

    Parameters:
    - tom_path (str): Path to the TOM matrix file (h5 format).
    - adata (AnnData): Annotated data object.
    - config (PlotConfig): Configuration object for plot control.
    - transcripts (List[str]): List of transcripts to keep in the network.
    - selected_module_colors (List[str]): List of module colors to display. If None, display all modules.
    - threshold (float): Threshold for TOM connection to consider an edge.
    - custom_filename (str): Custom filename for the plot.
    """

    TOM = load_subset_from_hdf5(
        filename=tom_path, rows=transcripts, cols=transcripts, threshold=threshold)
    module_info = adata.var

    if selected_module_colors:
        mod_genes = module_info[module_info['moduleColors'].isin(
            selected_module_colors)].index
    else:
        mod_genes = module_info.index

    mod_genes = [gene for gene in mod_genes if gene in TOM.index]

    G = nx.Graph()

    for i in range(len(mod_genes)):
        for j in range(i+1, len(mod_genes)):
            gene_i = mod_genes[i]
            gene_j = mod_genes[j]
            if TOM.loc[gene_i, gene_j] > threshold:
                G.add_edge(gene_i, gene_j, weight=TOM.loc[gene_i, gene_j])

    pos = nx.spring_layout(G, k=0.1, iterations=50)
    node_colors = [module_info.loc[node, 'moduleColors'] for node in G.nodes()]
    edge_weights = [G[u][v]['weight'] for u, v in G.edges()]

    plt.figure(figsize=(20, 20))
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=50)
    nx.draw_networkx_edges(G, pos, width=edge_weights)
    species = f"{adata.uns['species']}" if adata.uns['species'] is not None else ""
    if len(species.split()) > 1:
        plt.title(
            f"Co-Expression Network of $\it{{{species.split()[0]}}}$ $\it{{{species.split()[1]}}}$")
    else:
        plt.title(f"Co-Expression Network")

    full_name = f"co_expression_network" if not custom_filename else custom_filename

    if config.save_plots:
        plt.savefig(
            f"{config.output_path}/{full_name}.{config.file_format}", dpi=config.dpi)
        plt.savefig(f"{config.output_path}/{full_name}.png", dpi=config.dpi)

    if config.save_raw:
        raw_data_filepath = f"{config.output_path}/{full_name}.csv"
        raw_data = TOM.loc[list(G.nodes()), list(G.nodes())]
        raw_data.to_csv(raw_data_filepath)

    if config.show:
        plt.show()


def generate_shape_dict(traces: List[str], species_list_path: str = "../info/species") -> Dict[str, str]:
    """
    Generates a dictionary mapping each trace to unique Plotly marker symbols, ensuring consistent species-symbol mapping if a species file is provided.

    Parameters:
    - traces (List[str]): List of traces (e.g. modules).
    - species_list_path (str): Path to the file containing the list of species. If None, mapping will be random.

    Returns:
    - Dict: Dictionary mapping each module color to a unique marker symbol.
    """

    marker_symbols = [
        'circle', 'square', 'triangle-up', 'diamond', 'cross', 'x',
        'triangle-down', 'triangle-left', 'triangle-right', 'pentagon',
        'hexagon', 'hexagram', 'star', 'star-triangle-up', 'star-triangle-down',
        'star-square', 'star-diamond', 'triangle-ne', 'triangle-se', 'bowtie', 'hourglass',
        'triangle-sw', 'triangle-nw', 'hexagon2', 'diamond-tall', 'diamond-wide',
        'circle-open', 'square-open',
        'diamond-open', 'cross-dot', 'x-dot', 'triangle-up-dot', 'triangle-down-dot',
        'triangle-left-dot', 'triangle-right-dot', 'pentagon-dot', 'hexagon-dot',
        'star-dot', 'star-triangle-up-dot', 'star-triangle-down-dot', 'star-square-dot'
    ]

    default_species_list = [
        "Pteridophyllum racemosum",
        "Papaver somniferum",
        "Thalictrum thalictroides",
        "Aquilegia caerulea",
        "Staphisagria picta",
        "Capnoides sempervirens",
        "Nigella damascena",
        "Hydrastis canadensis",
        "Eschscholzia californica",
        "Epimedium grandiflorum"
    ]

    species_shape_dict = {}

    try:
        if species_list_path:
            with open(species_list_path, 'r') as file:
                species_list = [line.strip()
                                for line in file.readlines() if line.strip()]
        else:
            raise FileNotFoundError("No species file path provided.")
    except (FileNotFoundError, IOError) as e:
        print(f"Warning: {e}. Using default species list.")
        species_list = default_species_list

    # Ensure consistent mapping of species to symbols
    species_shape_dict = {
        species: marker_symbols[i % len(marker_symbols)]
        for i, species in enumerate(species_list)
    }

    unique_traces = list(set(traces))
    if species_shape_dict:
        shape_dict = {
            trace: species_shape_dict.get(
                trace, marker_symbols[i % len(marker_symbols)])
            for i, trace in enumerate(unique_traces)
        }
    else:
        shape_dict = {
            trace: marker_symbols[i % len(marker_symbols)]
            for i, trace in enumerate(unique_traces)
        }

    if len(unique_traces) > len(marker_symbols):
        print(
            f"Warning: More than {len(marker_symbols)} unique traces detected. Some traces may share the same marker symbol.")

    return shape_dict


def generate_cyto_shape_dict(traces: List[str], species_list_path: str = "../info/species") -> Dict[str, str]:
    """
    Generates a dictionary mapping each trace to Cytoscape.js-compatible shapes.

    Parameters:
    - traces (List[str]): List of traces (e.g. modules).
    - species_list_path (str): Path to the file containing the list of species. If None, mapping will be random.

    Returns:
    - Dict: Dictionary mapping each module color to a unique Cytoscape.js shape.
    """

    cyto_shapes = [
        'ellipse', 'triangle', 'rectangle', 'roundrectangle', 'cutrectangle',
        'bottomroundrectangle', 'barrel', 'diamond', 'pentagon', 'hexagon',
        'concavehexagon', 'heptagon', 'octagon', 'star', 'tag', 'vee'
    ]

    species_shape_dict = {}
    if species_list_path:
        # Read species from the file
        with open(species_list_path, 'r') as file:
            species_list = [line.strip() for line in file.readlines()]

        # Ensure consistent mapping of species to symbols
        species_shape_dict = {species: cyto_shapes[i % len(
            cyto_shapes)] for i, species in enumerate(species_list)}

    unique_traces = list(set(traces))
    if species_shape_dict:
        shape_dict = {trace: species_shape_dict.get(trace, cyto_shapes[i % len(
            cyto_shapes)]) for i, trace in enumerate(unique_traces)}
    else:
        shape_dict = {trace: cyto_shapes[i % len(
            cyto_shapes)] for i, trace in enumerate(unique_traces)}

    if len(unique_traces) > len(cyto_shapes):
        print(
            f"Warning: More than {len(cyto_shapes)} unique traces detected. Some traces may have the same shape.")

    return shape_dict


def plot_co_expression_network(tom: pd.DataFrame, adata: AnnData, config: PlotConfig,
                               selected_module_colors: List[str] = None, reduction_percentage: float = 0,
                               threshold: float = 0.5, colorscale: str = 'Viridis', custom_filename: str = None,
                               use_shapes: bool = False, use_colors: bool = False, title_suffix: str = None,
                               show_symbol_legend: bool = False, node_size: int = 10, show_edge_legend: bool = False,
                               cluster_name: str = None, additional_hover: str = None, width: int = 1200,
                               height: int = 800, identify_clusters: str = None, node_threshold_percent: float = 0.02,
                               highlight: List[str] = None, highlight_color: str = "magenta") -> None:
    """
    Plots the co-expression network for selected modules.

    Parameters:
    - tom (pd.DataFrame): The TOM matrix.
    - adata (AnnData): Annotated data object.
    - config (PlotConfig): Configuration object for plot control.
    - selected_module_colors (List[str]): List of module colors to display. If None, display all modules.
    - reduction_percentage (float): Percentage of genes to remove from the TOM matrix.
    - threshold (float): Threshold for TOM connection to consider an edge.
    - colorscale (str): The colorscale to use for the edge colors.
    - custom_filename (str): Custom filename for the plot.
    - use_shapes (bool): Use different shapes for different modules.
    - use_colors (bool): Use colors for modules.
    - title_suffix (str): Suffix to append to the plot title.
    - show_symbol_legend (bool): If True, display a legend for the node shapes.
    - node_size (int): Size of the nodes in the network.
    - show_edge_legend (bool): If True, display a legend for the edge colors.
    - cluster_name (str): Column name in adata.var for coloring nodes by cluster.
    - additional_hover (str): Column name from adata.var to include in hover text.
    - width (int): Width of the plot.
    - height (int): Height of the plot.
    - identify_clusters (str): If set, identify clusters in the network.
    - node_threshold_percent (float): Percentage of nodes required in a cluster to be assigned.
    - highlight (List[str]): List of nodes (e.g. transcripts) to highlight in the network.
    - highlight_color (str): Color to use for highlighting nodes.
    """

    module_info = adata.var

    if reduction_percentage > 0:
        tom = reduce_tom_matrix(tom, reduction_percentage)

    if selected_module_colors:
        mod_genes = module_info[module_info['moduleColors'].isin(
            selected_module_colors)].index
    else:
        mod_genes = module_info.index

    mod_genes = [gene for gene in mod_genes if gene in tom.index]

    G = nx.Graph()

    for i in range(len(mod_genes)):
        for j in range(i+1, len(mod_genes)):
            gene_i = mod_genes[i]
            gene_j = mod_genes[j]
            if tom.loc[gene_i, gene_j] > threshold:
                G.add_edge(gene_i, gene_j, weight=tom.loc[gene_i, gene_j])

    nodes_with_edges = set(G.nodes())

    cluster_map = None
    if identify_clusters and not cluster_name:
        cluster_map = identify_network_clusters(
            G, cluster_name=identify_clusters, node_threshold_percent=node_threshold_percent)
        cluster_name = identify_clusters
    elif identify_clusters and cluster_name:
        print("Warning: Both identify_clusters and cluster_name are set. Ignoring identify_clusters.")

    # First pass of spring layout
    # pos = nx.spring_layout(G, k=5, iterations=100)
    pos = nx.spring_layout(G, k=0.2, iterations=100)

    # Manually adjust positions of weakly connected nodes
    weakly_connected_nodes = [
        node for node in G.nodes() if G.degree(node) <= 2]
    for node in weakly_connected_nodes:
        neighbors = list(G.neighbors(node))
        if neighbors:
            neighbor_pos = np.mean([pos[neighbor]
                                   for neighbor in neighbors], axis=0)
            pos[node] = 0.5 * pos[node] + 0.8 * \
                neighbor_pos  # Adjust this factor as needed

    if len(G.edges or nodes_with_edges) == 0:
        fig = go.Figure()
        fig.update_layout(title="There is nothing to display",
                          xaxis=dict(showgrid=False, zeroline=False,
                                     showticklabels=False),
                          yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
        if config.show:
            fig.show()

        file_name = "co_expression_network" if not custom_filename else custom_filename

        if config.save_plots:
            if custom_filename:
                fig.write_image(
                    f"{config.output_path}/{custom_filename}.{config.file_format}", scale=config.dpi/96)
                fig.write_html(f"{config.output_path}/{custom_filename}.html")
            else:
                fig.write_image(
                    f"{config.output_path}/{file_name}.{config.file_format}", scale=config.dpi/96)
                fig.write_html(f"{config.output_path}/{file_name}.html")
        return

    edge_colors = [tom.loc[edge[0], edge[1]]
                   for edge in G.edges() if tom.loc[edge[0], edge[1]] > threshold]

    if edge_colors:
        cmin = min(edge_colors)
        cmax = max(edge_colors)

        if cmin == cmax:
            print(
                "Warning: All edge colors have the same value. Normalized edge colors will be set to 0.5.")
            normalized_edge_colors = [0.5] * len(edge_colors)
        else:
            normalized_edge_colors = [
                (color - cmin) / (cmax - cmin) for color in edge_colors]

        if any(np.isnan(color) for color in normalized_edge_colors):
            print("Error: Normalized edge colors contain NaN values.")
            return

        # Create edge traces with Viridis colors
        edge_traces = []
        for i, edge in enumerate(G.edges()):
            if tom.loc[edge[0], edge[1]] > threshold:
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                color = sample_colorscale(
                    colorscale, normalized_edge_colors[i])[0]
                edge_trace = go.Scatter(
                    x=[x0, x1, None],
                    y=[y0, y1, None],
                    line=dict(width=0.5, color=color),
                    hoverinfo='none',
                    mode='lines',
                    showlegend=False
                )
                edge_traces.append(edge_trace)
    else:
        edge_traces = []

    remaining_module_colors = list(
        module_info.loc[list(nodes_with_edges), 'moduleColors'])
    shape_dict = generate_shape_dict(remaining_module_colors)
    node_color_dict = {
        node: module_info.loc[node, 'moduleColors'] for node in nodes_with_edges}

    if cluster_map:
        cluster_strings = pd.Series(cluster_map).reindex(nodes_with_edges)

        def extract_cluster_id(cluster_string):
            try:
                return int(cluster_string.split(' ')[-1])
            except ValueError:
                return 0

        cluster_ids = cluster_strings.apply(
            extract_cluster_id).dropna().astype(int)
        cluster_colors = cluster_ids[cluster_ids != 0]

        # keep only nodes with cluster information
        nodes_with_cluster = cluster_colors.index
        unique_clusters = sorted(cluster_colors.unique())

        # Calculate average positions for clusters
        cluster_positions = {cluster: np.mean([pos[node] for node in nodes_with_cluster if cluster_colors[node] == cluster], axis=0)
                             for cluster in unique_clusters}

        # Apply an offset to separate clusters further apart
        offset_magnitude = 1.5  # Adjust this value to increase or decrease separation
        cluster_offsets = {cluster: np.array([np.random.uniform(-offset_magnitude, offset_magnitude),
                                              np.random.uniform(-offset_magnitude, offset_magnitude)])
                           for cluster in unique_clusters}

        for node in nodes_with_cluster:
            cluster = cluster_colors[node]
            pos[node] += cluster_offsets[cluster]

    # Update edge positions
    updated_edge_traces = []
    for i, edge in enumerate(G.edges()):
        if tom.loc[edge[0], edge[1]] > threshold:
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            color = sample_colorscale(colorscale, normalized_edge_colors[i])[0]
            edge_trace = go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                line=dict(width=0.5, color=color),
                hoverinfo='none',
                mode='lines',
                showlegend=False
            )
            updated_edge_traces.append(edge_trace)

    node_x = []
    node_y = []
    node_text = []
    node_colors = []
    node_symbols = []

    for node in nodes_with_edges:
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        hover_text = f"<b>Transcript:</b> {node}<br><b>Module:</b> {module_info.loc[node, 'moduleColors']}<br><b>Degree:</b> {len(list(G.neighbors(node)))}"
        if additional_hover and additional_hover in module_info.columns:
            hover_text += f"<br><b>{additional_hover}:</b> {module_info.loc[node, additional_hover]}"
        node_text.append(hover_text)

        if highlight and node in highlight:
            node_colors.append(highlight_color)  # Highlight nodes
        elif use_colors:
            node_colors.append(module_info.loc[node, 'moduleColors'])
        else:
            node_colors.append('black')

        if use_shapes:
            node_symbols.append(shape_dict.get(
                node_color_dict[node], 'circle'))
        else:
            node_symbols.append('circle')

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        text=node_text,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=False,
            size=node_size,
            color=node_colors,
            symbol=node_symbols,
        ),
        showlegend=False  # Ensure nodes themselves do not show up in legend
    )

    legend_shapes = []
    if show_symbol_legend and use_shapes:
        for color, shape in shape_dict.items():
            legend_shapes.append(go.Scatter(
                x=[None], y=[None],
                mode='markers',
                marker=dict(
                    size=node_size,
                    symbol=shape,
                    color='black'
                ),
                legendgroup=color,
                name=color,
                showlegend=True
            ))

    title = f"Co-Expression Network of <i>{adata.uns['species']}</i>"
    if title_suffix:
        title += f": {title_suffix}"

    if cluster_name:
        cluster_trace = go.Scatter(
            x=[cluster_positions[cluster][0] + cluster_offsets[cluster][0]
                for cluster in unique_clusters],
            y=[cluster_positions[cluster][1] + cluster_offsets[cluster][1]
                for cluster in unique_clusters],
            text=[str(cluster) for cluster in unique_clusters],
            mode='text',
            textposition='top center',
            textfont=dict(size=35, color='red', family='Arial', weight='bold'),
            hoverinfo='none',
            showlegend=False,
        )
        traces = updated_edge_traces + \
            [node_trace, cluster_trace] + legend_shapes
    else:
        traces = updated_edge_traces + [node_trace] + legend_shapes

    fig = go.Figure(data=traces,
                    layout=go.Layout(
                        title=dict(
                            text=title,
                            font=dict(size=20)
                        ),
                        showlegend=show_symbol_legend,
                        hovermode='closest',
                        xaxis=dict(showgrid=False, zeroline=False,
                                   showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )

    colorbar_trace = go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(
            colorscale='Viridis',
            cmin=cmin,
            cmax=cmax,
            colorbar=dict(
                title="TOM Similarity",
                thickness=20,
                len=0.5,
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.02,
                ticks="outside"
            )
        ),
        showlegend=False
    )

    if show_edge_legend:
        fig.add_trace(colorbar_trace)

    fig.update_layout(
        title={
            'x': 0.5,
            'xanchor': 'center'
        },
        width=width,
        height=height,
        margin=dict(b=20, l=20, r=20, t=40),
        legend=dict(
            x=-0.25,
            y=1,
            traceorder='normal',
            font=dict(size=15),
            itemsizing='constant',
            xanchor='left'
        ) if show_symbol_legend else None,
    )

    if config.show:
        fig.show()

    full_name = f"co_expression_network" if not custom_filename else custom_filename
    if config.save_plots:
        fig.write_image(f"{config.output_path}/{full_name}.{config.file_format}",
                        scale=config.dpi/96, width=1200, height=800)
        fig.write_html(f"{config.output_path}/{full_name}.html")

    if config.save_raw:
        raw_data_filepath = f"{config.output_path}/{full_name}.csv"
        raw_data = tom.loc[list(nodes_with_edges), list(nodes_with_edges)]
        raw_data.to_csv(raw_data_filepath)

    if identify_clusters:
        return cluster_map


#############
"""
ANNDATA PLOTS
Plotting functions for AnnData objects using calculated modules from WGCNA.
"""
#############


def plot_goea(adata: AnnData, config: PlotConfig, transcripts: Dict[str, List[str]], name: str = "PS",
              obo_path: str = "../go-basic.obo", go_column: str = "go_terms", top_percentage: int = 1,
              custom_filename: str = None) -> pd.DataFrame:
    """
    Performs and plots a GOEA using the top transcripts based on average expression.

    Parameters:
    - adata (AnnData): The AnnData object containing transcripts and their associated GO terms.
    - config (PlotConfig): Configuration object for plot control.
    - transcripts (Dict[str, List[str]]): Dictionary with species as keys and lists of transcripts.
    - obo_path (str): Path to the .obo file containing the GO ontology.
    - go_column (str): Column in adata.var containing comma separated GO terms for each transcript.
    - top_percentage (int): The top percentage of transcripts to be used for GOEA based on average expression.
    - custom_filename (str): Custom filename for the plot.

    Returns:
    - pd.DataFrame: A DataFrame containing significant GO terms and their details.
    """

    if top_percentage < 1:
        top_percentage = 1
        print("Top percentage set to 1 to ensure at least one transcript is selected.")

    go_dag = GODag(obo_path)

    # Get transcripts for the provided species
    transcripts = transcripts.get(adata.uns['species'], [])

    # Filter the adata.var to include only the provided transcript_ids
    filtered_var = adata.var.loc[transcripts]

    # Get average expression (mean_counts) for each transcript in transcript_ids
    avg_expression = filtered_var['mean_counts']

    # Determine the number of top transcripts
    num_top_transcripts = int(len(transcripts) * (top_percentage / 100))

    # Get top transcript IDs based on average expression
    top_transcripts = avg_expression.nlargest(num_top_transcripts).index
    top_transcript_ids = [
        t_id for t_id in transcripts if t_id in top_transcripts]

    # Prepare the GO background set
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
        methods=['fdr_bh']
    )

    # Perform GOEA
    goea_results = goea_obj.run_study(top_transcript_ids)

    # Filter significant results and extract additional information
    goea_results_sig = [
        (r.GO, r.name, r.p_fdr_bh,
         r.ratio_in_study[0]/r.ratio_in_study[1], r.ratio_in_study[0], go_dag[r.GO].depth)
        for r in goea_results if r.p_fdr_bh < 0.05
    ]

    if goea_results_sig:
        if config.save_plots:
            # Plot significant GO terms
            if custom_filename:
                plot_path = os.path.join(
                    config.output_path, f'{custom_filename}.{config.file_format}')
            else:
                plot_path = os.path.join(
                    config.output_path, f'goea_plot_{name}.{config.file_format}')

            significant_go_terms = [r[0] for r in goea_results_sig]
            plot_gos(plot_path, significant_go_terms, go_dag,
                     title=f'Significant GO Terms for {name}')

        # Create a DataFrame for significant results
        df_goea_results = pd.DataFrame(goea_results_sig, columns=[
                                       "GO_ID", "Term", "P-Value", "Fold Enrichment", "Num Genes", "Depth"])

        if config.save_raw:
            raw_data_filename = f"goea_results_{name}"
            raw_data_filepath = f"{config.output_path}/{raw_data_filename}.csv"
            df_goea_results.to_csv(raw_data_filepath, index=False)

        return df_goea_results
    else:
        print("No significant GO terms found.")
        return pd.DataFrame(columns=["GO_ID", "Term", "P-Value", "Fold Enrichment", "Num Genes", "Depth"])


def plot_highest_expr_genes(adata: AnnData, config: PlotConfig, n_top: int = 10, ortho: bool = False, custom_filename: str = None) -> None:
    """
    Wraps the scanpy.pl.highest_expr_genes function to add display and save functionality via a config object.

    Parameters:
    - adata (AnnData): The AnnData object containing the gene expression data.
    - config (PlotConfig): Configuration object for plot control.
    - n_top (int): Number of top expressing genes to display.
    - ortho (bool): If True, display orthogroups instead of transcripts.
    - custom_filename (str): Custom filename for the plot.
    """

    ax = sc.pl.highest_expr_genes(adata, n_top=n_top, show=False)
    fig = ax.get_figure()

    species = f"{adata.uns['species']}" if adata.uns['species'] is not None else ""
    identifier = "Orthogroups" if ortho else "Transcripts"

    if species != "":
        plt.title(
            f"Top {n_top} Highest Expressed {identifier} of $\it{{{species.split()[0]}}}$ $\it{{{species.split()[1]}}}$")
    else:
        plt.title(f"Top {n_top} Highest Expressed {identifier}")

    if config.show:
        plt.show()

    if config.save_plots:
        species = f"_{species}" if species != "" else ""
        if custom_filename:
            fig.savefig(
                f"{config.output_path}/{custom_filename}.{config.file_format}", dpi=config.dpi)
        else:
            fig.savefig(
                f"{config.output_path}/highest_expr_genes{species}.{config.file_format}", dpi=config.dpi)

    plt.close(fig)


def plot_violin(adata: AnnData, keys: List[str], config: PlotConfig, jitter: float = 0.4, groupby: str = None, custom_filename: str = None) -> None:
    """
    Wraps the scanpy.pl.violin function to add display and save functionality via a config object.

    Parameters:
    - adata (AnnData): The AnnData object containing the gene expression data.
    - keys (List[str]): Variables for which to plot the violin plots.
    - config (PlotConfig): Configuration object for plot control.
    - jitter (float): Amount of jitter to apply to the violin plot.
    - groupby (str): Key for grouping data points in the plot.
    - custom_filename (str): Custom filename for the plot.
    """

    # Generate the plot and capture the figure
    ax = sc.pl.violin(adata, keys, jitter=jitter, groupby=groupby, show=False)
    fig = ax.get_figure()

    if config.show:
        plt.show()

    if config.save_plots:
        if custom_filename:
            fig.savefig(
                f"{config.output_path}/{custom_filename}.{config.file_format}", dpi=config.dpi)
        else:
            fig.savefig(
                f"{config.output_path}/violin_plot_{'_'.join(keys)}.{config.file_format}", dpi=config.dpi)

    plt.close(fig)


def plot_rank_genes_groups_matrixplot(adata: AnnData, config: PlotConfig, n_genes: int = 4, groupby: str = "tissue", key: str = "dge_key",
                                      vmin: int = -5, vmax: int = 5, min_logfoldchange: int = 3, custom_filename: str = None) -> None:
    """
    Wraps the scanpy.pl.rank_genes_groups_matrixplot function to add display and save functionality via a config object.

    Parameters:
    - adata (AnnData): The AnnData object containing the gene expression data.
    - config (PlotConfig): Configuration object for plot control.
    - n_genes (int): Number of genes to display.
    - groupby (str): Column name to group by.
    - key (str): Key for differential gene expression results in adata.uns.
    - vmin (int): Minimum value for the color scale.
    - vmax (int): Maximum value for the color scale.
    - min_logfoldchange (int): Minimum log fold change to display.
    - custom_filename (str): Custom filename for the plot.
    """

    if groupby not in adata.obs.columns:
        groupby = "Combined_Trait"
        if groupby not in adata.obs.columns:
            raise ValueError(f"Column {groupby} not found in adata.obs.")

    ax_dict = sc.pl.rank_genes_groups_matrixplot(
        adata,
        n_genes=n_genes,
        groupby=groupby,
        key=key,
        values_to_plot="logfoldchanges",
        vmin=vmin,
        vmax=vmax,
        min_logfoldchange=min_logfoldchange,
        colorbar_title='log fold change',
        show=False
    )

    for ax in ax_dict.values():
        fig = ax.get_figure()

        if config.show:
            plt.show()

        if config.save_plots:
            if custom_filename:
                fig.savefig(
                    f"{config.output_path}/{custom_filename}.{config.file_format}", dpi=config.dpi)
            else:
                fig.savefig(
                    f"{config.output_path}/matrixplot_{groupby}_{key}.{config.file_format}", dpi=config.dpi)

        plt.close(fig)


def plot_correlation_network(df: pd.DataFrame, config: PlotConfig, threshold: float = 0.1, pvals: pd.DataFrame = None,
                             pval_threshold: float = 0.05, seed: int = None, color_map: str = "",
                             inter_species_only: bool = False, suffix: str = "", custom_filename: str = None,
                             positive: bool = True, use_shapes: bool = True, weight_label: str = "Edge Weights",
                             show_symbol_legend: bool = False, title: str = "Correlation Network",
                             width: int = 1200, height: int = 800, colorscale: str = 'Viridis',
                             show_edge_legend: bool = False) -> None:
    """
    Plots a correlation network based on a given dataframe. Optionally displays only inter-species edges.

    Parameters:
    - df (pd.DataFrame): The dataframe containing the correlation values.
    - config (PlotConfig): Configuration object for plot control.
    - threshold (float): The threshold value for including edges in the network.
    - pvals (pd.DataFrame): The dataframe containing the p-values.
    - pval_threshold (float): The threshold value for p-values to include edges in the network.
    - seed (int): The seed value for the random layout of the network.
    - color_map (str): If 'Modules', colors nodes based on their name prefixes; otherwise uses default coloring.
    - inter_species_only (bool): If True, only display edges between different species.
    - suffix (str): A suffix to append to the filename.
    - custom_filename (str): Custom filename for the plot.
    - positive (bool): If True, include only positive correlations; if False, include only negative correlations.
    - use_shapes (bool): If True, use different shapes for different species; otherwise use only one shape.
    - weight_label (str): The label for the edge weights.
    - show_symbol_legend (bool): If True, display a legend for the node shapes.
    - title (str): The title of the plot.
    - width (int): Width of the plot.
    - height (int): Height of the plot.
    - colorscale (str): The colorscale to use for the edge colors.
    - show_edge_legend (bool): If True, display a legend for the edge colors.
    """

    # Check if the DataFrame is symmetric and set diagonal to zero
    if df.equals(df.T):
        np.fill_diagonal(df.values, 0)

    if not use_shapes:
        show_symbol_legend = False

    G = nx.Graph()

    # Add edges between nodes with weights if they are above the threshold
    for source in df.index:
        # Get species name from the node label
        species_source = source.split('_')[-1]
        for target in df.columns:
            # Get species name from the node label
            species_target = target.split('_')[-1]
            weight = df.at[source, target]

            # Apply the correlation filter
            if positive and weight <= 0:
                continue
            if not positive and weight >= 0:
                continue

            # Apply the p-value filter if pvals is provided
            if pvals is not None:
                pval = pvals.at[source, target]
                if pval > pval_threshold:
                    continue

            if abs(weight) > threshold:
                if inter_species_only:
                    if species_source != species_target:
                        G.add_edge(source, target, weight=weight)
                else:
                    G.add_edge(source, target, weight=weight)

    # Remove isolates from the graph
    G.remove_nodes_from(list(nx.isolates(G)))

    if len(G.edges) == 0:
        print("No edges found above the threshold.")
        fig = go.Figure()
        fig.update_layout(title="There is nothing to display",
                          xaxis=dict(showgrid=False, zeroline=False,
                                     showticklabels=False),
                          yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
        if config.show:
            fig.show()

        file_name = f"correlation_network_{suffix}" if suffix else "correlation_network"

        if config.save_plots:
            if custom_filename:
                fig.write_image(
                    f"{config.output_path}/{custom_filename}.{config.file_format}", scale=config.dpi/96)
                fig.write_html(f"{config.output_path}/{custom_filename}.html")
            else:
                fig.write_image(
                    f"{config.output_path}/{file_name}.{config.file_format}", scale=config.dpi/96)
                fig.write_html(f"{config.output_path}/{file_name}.html")
        return

    # Handle edge weights
    weights = np.array([data['weight'] for _, _, data in G.edges(data=True)])
    if weights.size > 0:
        max_weight = weights.max()
        min_weight = weights.min()

        if min_weight == max_weight:
            scaled_weights = np.ones_like(
                weights) * 0.5  # Set to a constant value
        else:
            scaled_weights = (weights - min_weight) / (max_weight - min_weight)

        scaled_weights = np.nan_to_num(scaled_weights, nan=0.5)
        edge_colors = sample_colorscale(colorscale, scaled_weights)
    else:
        scaled_weights = np.array([])
        edge_colors = []

    pos = nx.spring_layout(G, k=1.5, iterations=200, seed=seed)

    edge_traces = []
    for idx, edge in enumerate(G.edges(data=True)):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        trace = go.Scatter(
            x=[x0, x1, None], y=[y0, y1, None],
            mode='lines',
            # Use color based on weight and set line width
            line=dict(width=1, color=edge_colors[idx]),
            hoverinfo='none',
            showlegend=False
        )
        edge_traces.append(trace)

    node_x = [pos[node][0] for node in G.nodes()]
    node_y = [pos[node][1] for node in G.nodes()]
    node_colors = []
    node_text = []
    node_symbols = []

    species_list = [node.split('_')[-1] for node in G.nodes()]
    if use_shapes:
        shape_dict = generate_shape_dict(species_list)
    else:
        shape_dict = {species: 'circle' for species in species_list}

    for node in G.nodes():
        species = node.split('_')[-1]
        symbol = shape_dict.get(species, 'circle')
        node_symbols.append(symbol)

        if color_map == 'Modules':
            parts = node.split('_')
            if len(parts) == 3:
                node_colors.append(parts[1])
            elif len(parts) == 2:
                node_colors.append(parts[0])
            else:
                parts = node.split('ME')
                if len(parts) == 2:
                    node_colors.append(parts[1])
        else:
            node_colors = 'black'

        node_degree = len(list(G.neighbors(node)))
        node_text.append(f"{node}<br>Degree: {node_degree}")

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        text=node_text,
        hoverinfo='text',
        marker=dict(
            symbol=node_symbols,
            showscale=False,
            size=15,
            color=node_colors,
            line_width=2),
        showlegend=False)

    edge_colorbar_trace = go.Scatter(
        x=[None],
        y=[None],
        mode='markers',
        marker=dict(
            colorscale=colorscale,
            cmin=min_weight,
            cmax=max_weight,
            colorbar=dict(
                title=weight_label,
                thickness=20,
                len=0.5,
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.02,
                ticks="outside"
            )
        ),
        hoverinfo='none',
        showlegend=False
    )

    legend_shapes = []
    if show_symbol_legend:
        if use_shapes:
            unique_species = set(species_list)
            for species in unique_species:
                shape = shape_dict.get(species, 'circle')
                legend_shapes.append(go.Scatter(
                    x=[None],
                    y=[None],
                    mode='markers',
                    marker=dict(
                        size=10,
                        symbol=shape,
                        color='black'
                    ),
                    legendgroup=species,
                    name=f"<i>{species}</i>",
                    showlegend=True
                ))

    fig = go.Figure(data=edge_traces + [node_trace] + legend_shapes,
                    layout=go.Layout(
                        title=dict(
                            text=title,
                            font=dict(size=20)
                        ),
                        showlegend=show_symbol_legend,
                        legend=dict(
                            x=-0.35,
                            y=1,
                            traceorder='normal',
                            font=dict(size=15),
                            itemsizing='constant',
                            xanchor='left'
                        ) if show_symbol_legend else None,
                        hovermode='closest',
                        margin=dict(b=20, l=20, r=20, t=40),
                        xaxis=dict(showgrid=False, zeroline=False,
                                   showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )

    if show_edge_legend:
        fig.add_trace(edge_colorbar_trace)

    fig.update_layout(
        title={
            'x': 0.5,
            'xanchor': 'center'
        },
        width=width,
        height=height
    )

    if config.show:
        fig.show()

    file_name = f"correlation_network_{suffix}" if suffix else "correlation_network"

    if config.save_plots:
        if custom_filename:
            fig.write_image(f"{config.output_path}/{custom_filename}.{config.file_format}",
                            scale=config.dpi/96, width=1200, height=800)
            fig.write_html(f"{config.output_path}/{custom_filename}.html")
        else:
            fig.write_image(f"{config.output_path}/{file_name}.{config.file_format}",
                            scale=config.dpi/96, width=1200, height=800)
            fig.write_html(f"{config.output_path}/{file_name}.html")

    if config.save_raw:
        raw_data_filename = custom_filename if custom_filename else f"correlation_network_raw_data_{suffix}"
        raw_data_filepath = f"{config.output_path}/{raw_data_filename}.csv"
        raw_data = df.loc[list(G.nodes()), list(G.nodes())]
        raw_data.to_csv(raw_data_filepath)


def add_legend_traces(min_weight: float, max_weight: float) -> List[go.Scatter]:
    """
    Adds traces for the colorbar legend based on the minimum and maximum edge weights.

    Parameters:
    - min_weight (float): The minimum edge weight.
    - max_weight (float): The maximum edge weight.

    Returns:
    - list: A list of Plotly traces for the colorbar legend.
    """

    legend_trace = []
    if max_weight > min_weight:
        # 8 different weight scales
        weight_scale = np.linspace(min_weight, max_weight, num=8)
        for weight in weight_scale:
            scaled_weight = (weight - min_weight) / \
                (max_weight - min_weight) * (1 - 0.1) + 0.1
            legend_trace.append(go.Scatter(
                x=[None], y=[None],
                mode='lines',
                line=dict(width=scaled_weight * 5, color='#888'),
                name=f'{weight:.2f}'
            ))
    else:
        standard_width = 5
        legend_trace.append(go.Scatter(
            x=[None], y=[None],
            mode='lines',
            line=dict(width=standard_width, color='#888'),
            name=f'{round(min_weight, 2)}'
        ))

    return legend_trace


def plot_overlap(overlap_matrix: pd.DataFrame, config: PlotConfig, column: str = "Tissues", custom_filename: str = None,
                 title: str = None, row_cluster: bool = True, col_cluster: bool = True, width: int = 1200, height: int = 800,
                 x_label: str = "Species", y_label: str = "Species") -> None:
    """
    Plots an interactive and responsive heatmap of the Jaccard similarity matrix of orthogroups across different species.
    The plot adjusts dynamically based on the screen size and remains interactive, allowing zoom and tooltip exploration.

    Parameters:
    - overlap_matrix (pd.DataFrame): A DataFrame containing Jaccard indices between each pair of combination.
    - config (PlotConfig): Configuration object for plot control.
    - column (str): The name to be used in the plot title.
    - custom_filename (str): Custom filename for the plot.
    - title (str): Title of the plot.
    - row_cluster (bool): Whether to cluster the rows.
    - col_cluster (bool): Whether to cluster the columns.
    - width (int): Width of the plot.
    - height (int): Height of the plot.
    - x_label (str): Label for the x-axis.
    - y_label (str): Label for the y-axis.
    """

    # Check if the overlap_matrix is empty or contains only NaN values
    if overlap_matrix.empty or overlap_matrix.isna().all().all():
        print("Warning: The provided overlap matrix is empty or contains only NaN values. Skipping plot.")
        return None

    # Perform clustering if enabled and the matrix is not empty
    if row_cluster:
        # Ensure there are enough rows for clustering
        if overlap_matrix.shape[0] > 1:
            row_linkage = linkage(overlap_matrix.values, method='ward')
            row_order = leaves_list(row_linkage)
            overlap_matrix = overlap_matrix.iloc[row_order, :]
        else:
            print("Warning: Not enough rows for row clustering. Skipping row clustering.")
            row_cluster = False

    if col_cluster:
        # Ensure there are enough columns for clustering
        if overlap_matrix.shape[1] > 1:
            col_linkage = linkage(overlap_matrix.T.values, method='ward')
            col_order = leaves_list(col_linkage)
            overlap_matrix = overlap_matrix.iloc[:, col_order]
        else:
            print(
                "Warning: Not enough columns for column clustering. Skipping column clustering.")
            col_cluster = False

    # Extract species from index and columns by splitting on '_'
    species_index = [idx.split('_')[-1] for idx in overlap_matrix.index]
    species_columns = [col.split('_')[-1] for col in overlap_matrix.columns]

    # Create a mask for intra-species comparisons (where species match across the index and columns)
    mask = np.array(species_index).reshape(-1, 1) == np.array(species_columns)

    # Create a copy of the matrix to apply mask
    masked_matrix = overlap_matrix.copy()
    masked_matrix[mask] = np.nan  # Mask intra-species comparisons as NaN

    # Determine the min and max of ALL values in the matrix, ignoring NaNs
    data_min_all = np.nanmin(
        overlap_matrix.values) if not overlap_matrix.isna().all().all() else 0
    data_max_all = np.nanmax(
        overlap_matrix.values) if not overlap_matrix.isna().all().all() else 1

    # Determine the min and max of the actual data range, ignoring both NaNs and the masked intra-species comparisons
    # Only consider values where the mask is False
    masked_values = overlap_matrix.values[~mask]
    if masked_values.size > 0:
        data_min = np.nanmin(masked_values)
        data_max = np.nanmax(masked_values)
    else:
        print("Warning: No valid values in the overlap matrix for plotting.")
        # Return empy plot if there are no valid values with "No valid values" title
        fig = go.Figure()
        fig.update_layout(title="No valid values")
        if config.show:
            fig.show()
        if config.save_plots:
            if custom_filename:
                fig.write_html(f"{config.output_path}/{custom_filename}.html")
            else:
                fig.write_html(f"{config.output_path}/overlap_{column}.html")
        if config.save_raw:
            raw_data_filename = custom_filename if custom_filename else f"overlap_{column}"
            raw_data_filepath = f"{config.output_path}/{raw_data_filename}.csv"
            overlap_matrix.to_csv(raw_data_filepath)

        return

    # Build the Plotly figure
    fig = go.Figure()

    # Add heatmap for actual data, scaling the color based on the actual data range
    fig.add_trace(go.Heatmap(
        z=masked_matrix.values,
        x=overlap_matrix.columns,
        y=overlap_matrix.index,
        colorscale='Viridis',
        zmin=data_min,
        zmax=data_max,
        # Colorbar adjustment
        colorbar=dict(
            title=dict(text="Jaccard Index", side="right"),
            # Set tick labels dynamically
            tickvals=[data_min, (data_min + data_max) / 2, data_max],
            # Tick formatting
            ticktext=[f'{data_min:.2f}',
                      f'{(data_min + data_max) / 2:.2f}', f'{data_max:.2f}'],
        ),
        hoverongaps=False,  # Don't show tooltips on masked areas
        texttemplate="%{z:.2f}",  # Display values inside the cells
        hoverinfo='x+y+z'
    ))

    # Add heatmap for masked intra-species comparisons (with a greyscale)
    # Use the original data for the masked areas to show the values
    fig.add_trace(go.Heatmap(
        # Only show values where species match
        z=np.where(mask, overlap_matrix.values, np.nan),
        x=overlap_matrix.columns,
        y=overlap_matrix.index,
        colorscale='Greys',
        zmin=data_min_all,
        zmax=data_max_all,  # Apply same scaling to maintain consistency
        hoverongaps=False,
        showscale=False,  # No colorbar for this greyscale heatmap
        texttemplate="%{z:.2f}",  # Show labels inside the cells
        hoverinfo='x+y+z'
    ))

    fig.update_layout(
        title=dict(
            text=title if title else f"Correlation of Orthogroups Between Different {column} and Species",
            font=dict(size=20),
            x=0.5,
            xanchor='center'
        ),
        xaxis_title=x_label,
        yaxis_title=y_label,
        autosize=True,
        width=width,
        height=height,
        margin=dict(l=150, r=10, b=150, t=80),
        hovermode='closest'
    )

    # Responsive resizing for axis labels and title
    fig.update_xaxes(tickangle=45, automargin=True)
    fig.update_yaxes(automargin=True)

    # Add zoom and pan
    fig.update_layout(
        dragmode='zoom',
        hoverlabel=dict(bgcolor="white", font_size=12)
    )

    if config.save_plots:
        if custom_filename:
            fig.write_html(f"{config.output_path}/{custom_filename}.html")
        else:
            fig.write_html(f"{config.output_path}/overlap_{column}.html")

    if config.show:
        fig.show()

    if config.save_raw:
        raw_data_filename = custom_filename if custom_filename else f"overlap_{column}"
        raw_data_filepath = f"{config.output_path}/{raw_data_filename}.csv"
        overlap_matrix.to_csv(raw_data_filepath)


def plot_stacked_results(results: Dict[str, pd.DataFrame], config: PlotConfig, custom_filename: str = None) -> Tuple[go.Figure, go.Figure]:
    """
    Visualizes stacked results for transcripts and GO terms per module from AnnData analysis.

    This function generates two stacked bar charts: one for the number of transcripts per module,
    and another for the number of unique GO terms per module. Each stack within a bar represents
    different categories from a specified .obs column if provided during the analysis. The hover
    text displays the category information.

    Parameters:
    - results (Dict[str, pd.DataFrame]): A dictionary of DataFrames from AnnData analysis. Each key represents a unique
      value from the .obs column or 'all_data' if obs_column is None. Each DataFrame contains
      'moduleColors' as the index, and two columns: 'transcripts_per_module' and 'unique_GO_terms_per_module'.
    - config (PlotConfig): Configuration object for plot control.
    - custom_filename (str): Custom filename for the plot.

    Returns:
    - Two Plotly figure objects: one for transcripts per module and another for unique GO terms per module.
    """

    # Diagram for transcripts per module
    fig_transcripts = go.Figure()
    # Diagram for unique GO terms per module
    fig_go_terms = go.Figure()

    for obs_value, df in results.items():
        # Adding trace for transcripts per module
        fig_transcripts.add_trace(go.Bar(
            name=obs_value,
            x=df['moduleColors'],
            y=df['transcripts_per_module'],
            hoverinfo='y+name'
        ))

        # Adding trace for unique GO terms per module
        fig_go_terms.add_trace(go.Bar(
            name=obs_value,
            x=df['moduleColors'],
            y=df['unique_GO_terms_per_module'],
            hoverinfo='y+name'
        ))

    # Update layouts for the figures
    fig_transcripts.update_layout(
        barmode='stack',
        title_text='Transcripts per Module',
        xaxis_title='Module Colors',
        yaxis_title='Number of Transcripts',
        hovermode='x unified'
    )

    fig_go_terms.update_layout(
        barmode='stack',
        title_text='Unique GO Terms per Module',
        xaxis_title='Module Colors',
        yaxis_title='Number of Unique GO Terms',
        hovermode='x unified'
    )

    if config.show:
        fig_transcripts.show()
        fig_go_terms.show()

    if config.save_plots:
        if custom_filename:
            fig_transcripts.write_image(
                f"{config.output_path}/{custom_filename}_transcripts.{config.file_format}", scale=config.dpi/96)
            fig_go_terms.write_image(
                f"{config.output_path}/{custom_filename}_go_terms.{config.file_format}", scale=config.dpi/96)
            fig_transcripts.write_html(
                f"{config.output_path}/{custom_filename}_transcripts.html")
            fig_go_terms.write_html(
                f"{config.output_path}/{custom_filename}_go_terms.html")
        else:
            fig_transcripts.write_image(
                f"{config.output_path}/transcripts.{config.file_format}", scale=config.dpi/96)
            fig_go_terms.write_image(
                f"{config.output_path}/go_terms.{config.file_format}", scale=config.dpi/96)
            fig_transcripts.write_html(
                f"{config.output_path}/transcripts.html")
            fig_go_terms.write_html(f"{config.output_path}/go_terms.html")

    return fig_transcripts, fig_go_terms


def plot_filtered_go_terms(df: pd.DataFrame, config: PlotConfig, p_value_threshold: float = 0.05, min_fold_enrichment: float = 0.2, min_depth: int = 2,
                           column: str = "Module", custom_filename: str = None) -> go.Figure:
    """
    Creates a bubble plot of Gene Ontology (GO) terms filtered by specified P-value threshold, minimum fold enrichment,
    and minimum depth. The size of each bubble represents the number of genes associated with the GO term, and each
    bubble is colored by the module to which the GO term belongs (if column is specified). Only GO terms with a P-value
    less than the specified threshold, a fold enrichment greater than the specified minimum, and a depth more than or
    equal to the specified minimum are included in the plot. The plot also displays the depth of each GO term when
    hovered over.

    Parameters:
    - df (pd.DataFrame): Containing the GO terms data. Expected columns are 'P-Value', 'Fold Enrichment', 'Num Genes', 'Module', 'GO_ID', and 'Depth'.
    - config (PlotConfig): Configuration object for plot control.
    - p_value_threshold (float): The maximum P-value for GO terms to be included in the plot.
    - min_fold_enrichment (float): The minimum fold enrichment for GO terms to be included in the plot.
    - min_depth (int): The minimum depth for GO terms to be included in the plot. If None, depth is not considered in the filtering.
    - column (str): The name of the column to use. If None, no column is used for coloring.
    - custom_filename (str): Custom filename for the plot.

    Returns:
    - go.Figure: A figure object containing the scatter plot.
    """

    if df.empty:
        print("No GO terms found.")
        fig = go.Figure()
        fig.update_layout(title="There is nothing to display",
                          xaxis=dict(showgrid=False, zeroline=False,
                                     showticklabels=False),
                          yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
        if config.show:
            fig.show()

        if config.save_plots:
            if custom_filename:
                fig.write_image(
                    f"{config.output_path}/{custom_filename}_{column if column else 'GO_Terms'}.{config.file_format}", scale=config.dpi/96)
                fig.write_html(
                    f"{config.output_path}/{custom_filename}_{column if column else 'GO_Terms'}.html")
            else:
                fig.write_image(
                    f"{config.output_path}/filtered_go_terms_{column if column else 'GO_Terms'}.{config.file_format}", scale=config.dpi/96)
                fig.write_html(
                    f"{config.output_path}/filtered_go_terms_{column if column else 'GO_Terms'}.html")
            return

    # Apply filters based on the specified thresholds and depth
    if min_depth is not None:
        filtered_df = df[(df['P-Value'] < p_value_threshold) & (df['Fold Enrichment']
                                                                >= min_fold_enrichment) & (df['Depth'] >= min_depth)].copy()
    else:
        filtered_df = df[(df['P-Value'] < p_value_threshold) &
                         (df['Fold Enrichment'] >= min_fold_enrichment)].copy()

    filtered_df['-log10(P-Value)'] = -np.log10(filtered_df['P-Value'])

    hover_data = {"Term": True,
                  "Depth": True,
                  "Num Genes": True,
                  }

    if column is not None:
        fig = px.scatter(filtered_df,
                         x="Fold Enrichment",
                         y="-log10(P-Value)",
                         size="Depth",
                         color=column,
                         hover_name="GO_ID",
                         hover_data=hover_data,
                         size_max=15,
                         title=f"GO Terms per {column}")

        fig.update_layout(xaxis_title='Fold Enrichment',
                          yaxis_title='-log10(P-Value)',
                          legend_title=column)

        first_trace = df[column].unique()[0]

        for trace in fig.data:
            trace.visible = 'legendonly' if trace.name != first_trace else True
    else:
        fig = px.scatter(filtered_df,
                         x="Fold Enrichment",
                         y="-log10(P-Value)",
                         size="Depth",
                         hover_name="GO_ID",
                         hover_data=hover_data,
                         size_max=15,
                         title="GO Terms")

        fig.update_layout(xaxis_title='Fold Enrichment',
                          yaxis_title='-log10(P-Value)')

    if config.show:
        fig.show()

    if config.save_plots:
        if custom_filename:
            fig.write_image(
                f"{config.output_path}/{custom_filename}_{column if column else 'GO_Terms'}.{config.file_format}", scale=config.dpi/96)
            fig.write_html(
                f"{config.output_path}/{custom_filename}_{column if column else 'GO_Terms'}.html")
        else:
            fig.write_image(
                f"{config.output_path}/filtered_go_terms_{column if column else 'GO_Terms'}.{config.file_format}", scale=config.dpi/96)
            fig.write_html(
                f"{config.output_path}/filtered_go_terms_{column if column else 'GO_Terms'}.html")

    if config.save_raw:
        raw_data_filename = custom_filename if custom_filename else f"filtered_go_terms_raw_data_{column if column else 'GO_Terms'}"
        raw_data_filepath = f"{config.output_path}/{raw_data_filename}.csv"
        filtered_df.to_csv(raw_data_filepath)

    return fig


def plot_umap(adata: AnnData, config: PlotConfig, obs_column: str = "tissue", custom_filename: str = None) -> go.Figure:
    """
    Creates an interactive UMAP scatter plot using Plotly, with points colored by the specified obs_column
    from an AnnData object. The specified column data is also shown when hovering over points.

    Parameters:
    - adata (AnnData): AnnData object. An AnnData object containing UMAP coordinates in adata.obsm['X_umap']
      and categorical data in adata.obs for coloring and hover information.
    - config (PlotConfig): Configuration object for plot control.
    - obs_column (str): The name of the column in adata.obs to use for coloring and hover data.
    - custom_filename (str): Custom filename for the plot.

    Returns:
    - go.Figure: A figure object containing the scatter plot.
    """

    if obs_column not in adata.obs.columns:
        obs_column = "Combined_Trait"
        if obs_column not in adata.obs.columns:
            raise ValueError(f"Column {obs_column} not found in adata.obs.")

    custom_colors = [
        '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896',
        '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7',
        '#bcbd22', '#dbdb8d', '#17becf', '#9edae5', '#393b79', '#5254a3', '#6b6ecf', '#9c9ede',
        '#637939', '#8ca252', '#b5cf6b', '#cedb9c', '#8c6d31', '#bd9e39'
    ]

    c_obs = obs_column.capitalize()

    umap_data = pd.DataFrame(
        adata.obsm['X_umap'], columns=['UMAP-1', 'UMAP-2'])
    umap_data[c_obs] = adata.obs[obs_column].values
    umap_data['Sample'] = adata.obs.index.values

    fig = px.scatter(umap_data, x='UMAP-1', y='UMAP-2',
                     color=c_obs,
                     color_discrete_map={value: custom_colors[i % 30] for i, value in enumerate(
                         adata.obs[obs_column].unique())},
                     hover_name='Sample',
                     hover_data={c_obs: True})

    fig.update_layout(
        legend_title_text=c_obs,
    )

    fig.update_traces(marker=dict(size=10))

    if config.show:
        fig.show()

    if config.save_plots:
        if custom_filename:
            fig.write_image(
                f"{config.output_path}/{custom_filename}.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/{custom_filename}.html")
        else:
            fig.write_image(
                f"{config.output_path}/umap_plot.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/umap_plot.html")

    return fig


def plot_pca_overview(adata: AnnData, config: PlotConfig, color: str = "tissue", size: int = 50, custom_filename: str = None) -> None:
    """
    Plots PCA overview with Scanpy, then manually manages plot saving for each generated figure.

    Parameters:
    - adata (AnnData): The AnnData object containing PCA results.
    - config (PlotConfig): Configuration object for plot control.
    - color (str): The observation key to color by.
    - size (int): The size of the points in the plot.
    - custom_filename (str): Custom filename for the plot.
    """

    if color not in adata.obs.columns:
        color = "Combined_Trait"
        if color not in adata.obs.columns:
            raise ValueError(f"Column {color} not found in adata.obs.")

    # Generate the PCA overview plots
    sc.pl.pca_overview(adata, color=color, show=False, size=size)

    # Loop through all open figures and save them
    if config.save_plots:
        for i, fig in enumerate(plt.get_fignums()):
            fig_obj = plt.figure(fig)
            if custom_filename:
                fig_obj.savefig(
                    f"{config.output_path}/{custom_filename}_fig{i}.{config.file_format}", dpi=config.dpi)
            else:
                fig_obj.savefig(
                    f"{config.output_path}/pca_overview_{color}_fig{i}.{config.file_format}", dpi=config.dpi)
            plt.close(fig)

    if config.show:
        sc.pl.pca_overview(adata, color=color, show=True, size=size)


def get_and_plot_explained_variance(adata: AnnData, config: PlotConfig, target_variance: float = 0.98, custom_filename: str = None) -> Tuple[int, go.Figure]:
    """
    Plots the explained variance per principal component (PC) and the cumulative explained variance using Plotly.

    This function takes an AnnData object containing PCA results and plots two lines:
    1. The explained variance by each PC.
    2. The cumulative explained variance to understand how many PCs are required to explain a certain percentage of the variance.

    A vertical line is also plotted to indicate the number of PCs required to reach a target of 90% explained variance.

    Parameters:
    - adata (AnnData): Must contain PCA results in adata.uns['pca']['variance_ratio'].
    - config (PlotConfig): Configuration object for plot control.
    - target_variance (float): The target variance to reach.
    - custom_filename (str): Custom filename for the plot.

    Returns:
    - Tuple[int, go.Figure]: The number of PCs required to reach the target variance and a figure object containing the bar plot.
    """

    # Extract the explained variance ratio per PC
    explained_variance_ratio = adata.uns['pca']['variance_ratio']

    # Calculate the cumulative explained variance
    cumulative_explained_variance = np.cumsum(explained_variance_ratio)

    # Calculate the number of PCs required to reach the target variance
    num_pcs_for_target_variance = next(i for i, cum_var in enumerate(
        cumulative_explained_variance, 1) if cum_var >= target_variance)

    # Create a plot with both explained and cumulative explained variance
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=list(range(1, len(explained_variance_ratio) + 1)),
        y=explained_variance_ratio,
        mode='lines+markers',
        name='Explained Variance per PC',
        hoverinfo='text',
        text=[f'<b>PC:</b> {i}<br><b>Explained Variance:</b> {var:.2%}' for i,
              var in enumerate(explained_variance_ratio, 1)]
    ))
    fig.add_trace(go.Scatter(
        x=list(range(1, len(cumulative_explained_variance) + 1)),
        y=cumulative_explained_variance,
        mode='lines+markers',
        name='Cumulative Explained Variance',
        hoverinfo='text',
        text=[f'<b>PC:</b> {i}<br><b>Cumulative Explained Variance:</b> {cum_var:.2%}' for i,
              cum_var in enumerate(cumulative_explained_variance, 1)]
    ))

    # Add a vertical line for the target variance
    fig.add_shape(type="line", x0=num_pcs_for_target_variance, x1=num_pcs_for_target_variance, y0=0, y1=target_variance,
                  line=dict(color="Red", width=2, dash="dash"))
    fig.add_annotation(x=num_pcs_for_target_variance, y=target_variance, text=f"PC {num_pcs_for_target_variance} reaches {target_variance * 100:.0f}% variance",
                       showarrow=True, arrowhead=1, ax=20)

    fig.update_layout(
        title='PCA - Explained Variance',
        xaxis_title='Number of Principal Components (PCs)',
        yaxis_title='Explained Variance',
        legend_title="Legend",
        showlegend=True
    )

    fig.update_xaxes(range=[0, len(explained_variance_ratio) + 1])
    fig.update_yaxes(range=[0, 1])

    if config.show:
        fig.show()

    if config.save_plots:
        if custom_filename:
            fig.write_image(
                f"{config.output_path}/{custom_filename}.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/{custom_filename}.html")
        else:
            fig.write_image(
                f"{config.output_path}/explained_variance.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/explained_variance.html")

    return num_pcs_for_target_variance, fig


def plot_module_gene_proportions(adata: AnnData, config: PlotConfig, custom_filename: str = None) -> go.Figure:
    """
    Create a sorted bar plot showing the proportion of genes in each module using Plotly.

    Parameters:
    - adata (AnnData): The annotated data matrix.
    - config (PlotConfig): Configuration object for plot control.
    - custom_filename (str): Custom filename for the plot.

    Returns:
    - go.Figure: A figure object containing the sorted bar plot.
    """

    total_genes = adata.var.shape[0]
    module_colors = adata.var['moduleColors'].unique()
    proportions = []

    for module in module_colors:
        module_genes = adata.var[adata.var['moduleColors'] == module].shape[0]
        proportion = round((module_genes / total_genes) * 100, 2)
        proportions.append((module, proportion))

    # Sort modules by proportion, descending
    sorted_modules = sorted(proportions, key=lambda x: x[1], reverse=True)
    sorted_module_colors = [module[0] for module in sorted_modules]
    sorted_proportions = [module[1] for module in sorted_modules]

    # Create sorted bar plot
    fig = go.Figure(data=[go.Bar(x=sorted_module_colors,
                    y=sorted_proportions, marker_color=sorted_module_colors)])
    fig.update_layout(
        title='Proportion of Transcripts in Each Module',
        xaxis_title='Module',
        yaxis_title='Percentage of Total Transcripts (%)',
        xaxis_tickangle=-45
    )

    if config.show:
        fig.show()

    if config.save_plots:
        if custom_filename:
            fig.write_image(
                f"{config.output_path}/{custom_filename}.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/{custom_filename}.html")
        else:
            fig.write_image(
                f"{config.output_path}/module_gene_proportions.{config.file_format}", scale=config.dpi/96)
            fig.write_html(
                f"{config.output_path}/module_gene_proportions.html")

    return fig


def plot_transcripts_and_go(data: Dict[str, pd.DataFrame], config: PlotConfig, custom_filename: str = None) -> None:
    """
    Create separate line charts for transcripts and GO terms for various tissue stages.
    Only the first tissue is displayed by default, others are available in the legend.

    Parameters:
    - data (Dict[str, pd.DataFrame]): A dictionary of DataFrames, where each key is a tissue stage and each value is a DataFrame containing
      columns 'moduleColors', 'transcripts_per_module', and 'unique_GO_terms_per_module'.
    - config (PlotConfig): Configuration object for plot control.
    - custom_filename (str): Custom filename for the plot.
    """

    fig_transcripts = go.Figure()
    fig_go_terms = go.Figure()

    # Assign colors to each line, one color per tissue stage
    colors = px.colors.qualitative.Plotly

    first = True  # Flag to track the first trace
    for i, (stage, df) in enumerate(data.items()):
        # Cyclically select a color from the palette
        color = colors[i % len(colors)]
        # Only the first trace is visible initially
        visible = True if first else 'legendonly'

        # Add lines for transcripts
        fig_transcripts.add_trace(go.Scatter(
            x=df['moduleColors'], y=df['transcripts_per_module'],
            mode='lines+markers',
            name=stage,
            line=dict(color=color),
            visible=visible
        ))

        # Add lines for GO terms
        fig_go_terms.add_trace(go.Scatter(
            x=df['moduleColors'], y=df['unique_GO_terms_per_module'],
            mode='lines+markers',
            name=stage,
            line=dict(color=color),
            visible=visible
        ))

        first = False

    fig_transcripts.update_layout(
        title='Transcripts per Module and Tissue Stage',
        xaxis=dict(title='Module Colors', tickangle=-45),
        yaxis=dict(title='Number of Transcripts'),
        legend=dict(title='Tissue Stages')
    )

    fig_go_terms.update_layout(
        title='GO Terms per Module and Tissue Stage',
        xaxis=dict(title='Module Colors', tickangle=-45),
        yaxis=dict(title='Number of GO Terms'),
        legend=dict(title='Tissue Stages')
    )

    if config.show:
        fig_transcripts.show()
        fig_go_terms.show()

    if config.save_plots:
        if custom_filename:
            fig_transcripts.write_image(
                f"{config.output_path}/{custom_filename}_transcripts.{config.file_format}", scale=config.dpi/96)
            fig_go_terms.write_image(
                f"{config.output_path}/{custom_filename}_go_terms.{config.file_format}", scale=config.dpi/96)
            fig_transcripts.write_html(
                f"{config.output_path}/{custom_filename}_transcripts.html")
            fig_go_terms.write_html(
                f"{config.output_path}/{custom_filename}_go_terms.html")
        else:
            fig_transcripts.write_image(
                f"{config.output_path}/transcripts.{config.file_format}", scale=config.dpi/96)
            fig_go_terms.write_image(
                f"{config.output_path}/go_terms.{config.file_format}", scale=config.dpi/96)
            fig_transcripts.write_html(
                f"{config.output_path}/transcripts.html")
            fig_go_terms.write_html(f"{config.output_path}/go_terms.html")


def plot_hub_connectivity(adata: AnnData, config: PlotConfig, top_n: int = 1, primary_y_column: str = 'connectivity', secondary_y_column: str = 'mean_counts',
                          primary_y_label: str = 'Connectivity', secondary_y_label: str = 'Mean Counts',
                          hover_columns: List[str] = ["ortho_ID", "connectivity", "mean_counts", "total_counts"], custom_filename: str = None) -> go.Figure:
    """
    Plots a dual y-axis bar plot for specified columns of gene data with customizable hover information.
    Accepts an AnnData object directly, ensuring that additional information from adata.var is included.

    Parameters:
    - adata (AnnData): AnnData object containing the bulk RNAseq data.
    - config (PlotConfig): Configuration object for plot control.
    - top_n (int): The number of top hub genes to include in the plot for each module.
    - primary_y_column (str): The name of the column for the primary y-axis.
    - secondary_y_column (str): The name of the column for the secondary y-axis.
    - primary_y_label (str): Label for the primary y-axis.
    - secondary_y_label (str): Label for the secondary y-axis.
    - hover_columns (List[str]): List of column names to include in the hover tooltips.
    - custom_filename (str): Custom filename for the plot.

    Returns:
    - go.Figure: A figure object containing the bar plot.
    """

    # Ensure secondary_y_column is included in hover data
    if secondary_y_column not in hover_columns:
        hover_columns.append(secondary_y_column)

    all_data = []
    for module, df in adata.uns['hub_genes'].items():
        df = df.nlargest(top_n, primary_y_column)
        df['moduleColor'] = module
        df['Transcript'] = df.index
        if 'gene_id' in df.columns:
            df = df.set_index('gene_id')
        needed_var_cols = [
            col for col in hover_columns if col not in df.columns and col in adata.var.columns]
        if needed_var_cols:
            df = df.join(adata.var[needed_var_cols], how='left')
        # Reset index to maintain gene_id if needed
        all_data.append(df.reset_index())

    combined_df = pd.concat(all_data).sort_values(
        by=primary_y_column, ascending=False)

    # Create hover text with rounded values and include the original index
    hover_text = combined_df.apply(lambda row: '<br>'.join([
        f'<b>Transcript</b>: {row["Transcript"]}',
        f'<b>Connectivity</b>: {round(float(row["connectivity"]), 2)}',
        f'<b>Mean counts</b>: {round(float(row["mean_counts"]), 2)}',
        f'<b>Total counts</b>: {round(float(row["total_counts"]), 2)}'] +
        [f'<b>{col}</b>: {row[col]}' for col in hover_columns if col not in [
            "connectivity", "mean_counts", "total_counts", "Transcript"]]
    ), axis=1)

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    x_positions = list(range(len(combined_df)))
    offset = 0.2
    bar_width = 0.35

    # Add primary y-axis bar plot
    fig.add_trace(
        go.Bar(x=[x - offset for x in x_positions], y=combined_df[primary_y_column],
               name=primary_y_label, marker=dict(color='DarkSlateGrey'), width=bar_width,
               hovertext=hover_text, hoverinfo='text'),
        secondary_y=False
    )

    # Add secondary y-axis bar plot
    fig.add_trace(
        go.Bar(x=[x + offset for x in x_positions], y=combined_df[secondary_y_column],
               name=secondary_y_label, marker=dict(color='Crimson'), width=bar_width,
               hovertext=hover_text, hoverinfo='text'),
        secondary_y=True
    )

    fig.update_layout(
        title=f'{primary_y_label} and {secondary_y_label} by Module',
        xaxis_title='Module Colors',
        xaxis=dict(
            tickmode='array',
            tickvals=x_positions,
            ticktext=combined_df['moduleColor']
        ),
        xaxis_tickangle=-45,
        barmode='group'
    )
    fig.update_yaxes(title_text=f'{primary_y_label}', secondary_y=False)
    fig.update_yaxes(title_text=f'{secondary_y_label}', secondary_y=True)

    if config.show:
        fig.show()

    if config.save_plots:
        if custom_filename:
            fig.write_image(
                f"{config.output_path}/{custom_filename}.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/{custom_filename}.html")
        else:
            fig.write_image(
                f"{config.output_path}/hub_connectivity.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/hub_connectivity.html")

    return fig


def plot_orthos_and_transcripts_per_module(adata: AnnData, config: PlotConfig, ortho_id_col: str = 'ortho_ID', module_color_col: str = 'moduleColors', custom_filename: str = None) -> go.Figure:
    """
    Creates a grouped bar plot showing the number of unique ortholog IDs and total number of transcripts per module color
    with detailed hover text.

    Parameters:
    adata (AnnData): The anndata object containing ortholog IDs and module colors.
    ortho_id_col (str): Column name for ortholog IDs. Default is 'ortho_ID'.
    module_color_col (str): Column name for module colors. Default is 'moduleColors'.
    config (PlotConfig): Configuration object for plot control.
    custom_filename (str): Custom filename for the plot.

    Returns:
    - Figure: A figure object containing the bar plot.
    """

    df_filtered = adata.var.dropna(subset=[ortho_id_col])
    ortho_counts = df_filtered.groupby(module_color_col)[
        ortho_id_col].nunique()
    module_color_counts = adata.var[module_color_col].value_counts()

    result_df = pd.DataFrame({
        'Number of Orthogroups': ortho_counts,
        'Number of Transcripts': module_color_counts
    }).fillna(0)

    result_df.sort_values('Number of Transcripts',
                          ascending=False, inplace=True)

    fig = make_subplots(specs=[[{"secondary_y": True}]])

    x_positions = list(range(len(result_df)))
    offset = 0.2  # offset for grouped bars
    bar_width = 0.35  # Smaller bar width for clearer separation

    # Generate hover text for each bar
    hover_text_orthos = ['<b>Module Color:</b> {}<br><b>Number of Orthogroups:</b> {}'.format(module, count)
                         for module, count in zip(result_df.index, result_df['Number of Orthogroups'])]
    hover_text_transcripts = ['<b>Module Color:</b> {}<br><b>Number of Transcripts:</b> {}'.format(module, count)
                              for module, count in zip(result_df.index, result_df['Number of Transcripts'])]

    # Add the Number of Orthogroups bar
    fig.add_trace(
        go.Bar(x=[x - offset for x in x_positions], y=result_df['Number of Orthogroups'],
               name='Number of Orthogroups', marker=dict(color='DarkSlateGrey'), width=bar_width,
               hoverinfo='text', hovertext=hover_text_orthos),
        secondary_y=False
    )

    # Add the Number of Transcripts bar
    fig.add_trace(
        go.Bar(x=[x + offset for x in x_positions], y=result_df['Number of Transcripts'],
               name='Number of Transcripts', marker=dict(color='Crimson'), width=bar_width,
               hoverinfo='text', hovertext=hover_text_transcripts),
        secondary_y=True
    )

    fig.update_layout(
        title='Orthology Groups and Transcripts by Module',
        xaxis_title='Module Colors',
        xaxis=dict(
            tickmode='array',
            tickvals=x_positions,
            ticktext=result_df.index
        ),
        xaxis_tickangle=-45,
        barmode='group'
    )

    fig.update_yaxes(title_text='Number of Orthogroups',
                     secondary_y=False)
    fig.update_yaxes(title_text='Number of Transcripts', secondary_y=True)

    if config.show:
        fig.show()

    if config.save_plots:
        if custom_filename:
            fig.write_image(
                f"{config.output_path}/{custom_filename}.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/{custom_filename}.html")
        else:
            fig.write_image(
                f"{config.output_path}/orthos_transcripts.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/orthos_transcripts.html")

    return fig


def plot_eigengenes_correlations(adata: AnnData, config: PlotConfig, custom_filename: str = None) -> None:
    """
    Calculate the correlation matrix of module eigen genes and plot a heatmap, 
    with the diagonal greyed out and the color scale adjusted to non-diagonal values.

    Parameters:
    - adata: The AnnData object containing the module eigengenes.
    - config: Configuration object for plot control.
    - custom_filename: Custom filename for the plot.
    """

    correlation_matrix = adata.obsm['eigengenes'].corr()

    mask = np.zeros_like(correlation_matrix, dtype=bool)
    np.fill_diagonal(mask, True)

    fig, ax = plt.subplots(figsize=(12, 10))
    max_scale_value = np.max(np.abs(correlation_matrix.values[~mask]))

    sns.heatmap(correlation_matrix, ax=ax, mask=mask, annot=False, cmap='coolwarm', fmt=".2f",
                linewidths=.5, center=0, vmax=max_scale_value, vmin=-max_scale_value,
                cbar_kws={'label': 'Correlation Coefficient'},
                linecolor='black', square=True)

    ax.set_title('Correlation between Module Eigengenes')

    for i in range(len(mask)):
        for j in range(len(mask)):
            if mask[i, j]:
                ax.add_patch(plt.Rectangle(
                    (j, i), 1, 1, fill=True, color='black'))

    if config.show:
        plt.show()

    if config.save_plots:
        if custom_filename:
            file_path = f"{config.output_path}/{custom_filename}.{config.file_format}"
            fig.savefig(file_path, dpi=config.dpi)
        else:
            file_path = f"{config.output_path}/eigengenes_correlations.{config.file_format}"
            fig.savefig(file_path, dpi=config.dpi)

    if config.save_raw:
        raw_data_filename = custom_filename if custom_filename else "eigengenes_correlation_matrix"
        raw_data_filepath = f"{config.output_path}/{raw_data_filename}.csv"
        correlation_matrix.to_csv(raw_data_filepath)

    plt.close(fig)


###########
"""
ORTHO PLOTS
Plotting functions for orthology data. 
Most functions are designed to work with the ortho_df DataFrame which can be generated using "get_ortho_count" function from utils module.
"""
###########


def plot_histogram(df: pd.DataFrame, config: PlotConfig, column: str = 'Count', bins: int = 40, custom_filename: str = None) -> None:
    """
    Plots a histogram of the specified column in the DataFrame.

    Parameters:
    - df (pd.DataFrame): The DataFrame to plot.
    - config (PlotConfig): The configuration object for plot control.
    - column (str): The column name to plot the histogram for.
    - bins (int): The number of bins in the histogram.
    - custom_filename (str): The custom filename for the plot.
    """

    plt.figure(figsize=(10, 6))
    plt.hist(df[column], bins=bins, color='skyblue', edgecolor='black')
    plt.title(f'{column} of Orthogroups')
    plt.xlabel(column)
    plt.ylabel('Frequency')

    if config.save_plots:
        if custom_filename:
            plt.savefig(
                f"{config.output_path}/{custom_filename}.{config.file_format}", dpi=config.dpi)
        else:
            plt.savefig(
                f"{config.output_path}/histogram_{column}.{config.file_format}", dpi=config.dpi)

    if config.show:
        plt.show()

    plt.close()


def plot_venn(df: pd.DataFrame, config: PlotConfig, species_list: List[str], custom_filename: str = None, italic_labels: bool = True, colors: List[str] = ['skyblue', 'orange', 'green'], alpha: float = 0.7) -> None:
    """
    Plots a Venn diagram for the given species list from the ortho_df DataFrame.

    Parameters:
    - df (DataFrame): DataFrame containing ortho data.
    - config (PlotConfig): Configuration object for plot control.
    - species_list (List[str]): List of three species to include in the Venn diagram.
    - custom_filename (str): Custom filename for the plot.
    - italic_labels (bool): Whether to italicize the labels in the Venn diagram.
    - colors (List[str]): List of colors for the Venn diagram.
    - alpha (float): Transparency of the Venn diagram patches.
    """

    if len(species_list) != 3:
        print("Venn diagram can only be plotted for exactly three species.")
        return

    sets = [set(df[df[species] > 0].index) for species in species_list]

    plt.figure(figsize=(10, 8))

    if colors:
        v = venn3(sets, set_labels=species_list, set_colors=colors)
        colors = [to_rgba(color, alpha=alpha) for color in colors]
    else:
        v = venn3(sets, set_labels=species_list)

    plt.title('Venn Diagram of Orthogroups', fontsize=20)

    for label in v.set_labels:
        if label:
            label.set_fontsize(15)
            if italic_labels:
                label.set_fontstyle('italic')

    if config.save_plots:
        if custom_filename:
            plt.savefig(
                f"{config.output_path}/{custom_filename}.{config.file_format}", dpi=config.dpi)
        else:
            plt.savefig(
                f"{config.output_path}/venn_diagram.{config.file_format}", dpi=config.dpi)

    if config.show:
        plt.show()

    plt.close()


def plot_upset(df: pd.DataFrame, config: PlotConfig, species_list: List[str], custom_filename: str = None, italic_labels: bool = True) -> None:
    """
    Plots an UpSet plot for the given species list using the prepared DataFrame.

    Parameters:
    - df (pd.DataFrame): DataFrame containing ortho data in the binary format expected by UpSetPlot.
    - config (PlotConfig): Configuration object for plot control.
    - species_list (List[str]): List of species to include in the UpSet plot.
    - custom_filename (str): Custom filename for the plot.
    - italic_labels (bool): Whether to italicize the species names in the category labels of the UpSet plot.
    """

    df_filtered = df[species_list]
    binary_df = (df_filtered > 0).astype(int)
    binary_df = binary_df.set_index(species_list)

    upset = UpSet(binary_df, subset_size='count', sort_by='cardinality')
    plot = upset.plot()
    plt.title('UpSet Plot of Orthogroups', fontsize=20)
    plt.ylabel('Number of Orthogroups', fontsize=12)

    if italic_labels:
        for label in plot['matrix'].get_yticklabels():
            label.set_fontstyle('italic')

    if config.save_plots:
        if custom_filename:
            plt.savefig(
                f"{config.output_path}/{custom_filename}.{config.file_format}", dpi=config.dpi)
        else:
            plt.savefig(
                f"{config.output_path}/upset_plot.{config.file_format}", dpi=config.dpi)

    if config.show:
        plt.show()

    plt.close()


def plot_pca(df: pd.DataFrame, config: PlotConfig, custom_filename: str) -> None:
    """
    Plots a PCA of the ortho_df DataFrame.

    Parameters:
    - df (pd.DataFrame): DataFrame to plot.
    - config (PlotConfig): Configuration object for plot control.
    - custom_filename (str): Custom filename for the plot.
    """

    data = df.drop(columns='Count').T
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(data)

    principal_df = pd.DataFrame(
        data=principal_components, columns=['PC1', 'PC2'])
    principal_df['Species'] = data.index

    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=principal_df, x='PC1', y='PC2',
                    hue='Species', palette='tab10')
    plt.title('PCA of Ortho Counts by Species')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.legend(title='Species', bbox_to_anchor=(1.05, 1), loc='upper left')

    if config.save_plots:
        if custom_filename:
            plt.savefig(
                f"{config.output_path}/{custom_filename}.{config.file_format}", dpi=config.dpi)
        else:
            plt.savefig(
                f"{config.output_path}/pca_plot.{config.file_format}", dpi=config.dpi)

    if config.show:
        plt.show()

    plt.close()


def plot_top_ortho_groups(df: pd.DataFrame, config: PlotConfig, top_n: int = 20, custom_filename: str = None, italic_labels: bool = True) -> None:
    """
    Plots a subset of the DataFrame containing only the top N ortho groups based on 'Count'.

    Parameters:
    - df (pd.DataFrame): DataFrame to plot.
    - config (PlotConfig): Configuration object for plot control.
    - top_n (int): Number of top records to plot.
    - custom_filename (str): Custom filename for the plot.
    - italic_labels (bool): Whether to italicize the ortho group labels in the plot.
    """

    # Sort and select the top N rows
    top_df = df.sort_values(by='Count', ascending=False).head(top_n)

    fig = go.Figure()

    for column in top_df.columns.drop('Count'):
        fig.add_trace(go.Bar(
            x=top_df.index,
            y=top_df[column],
            name=column
        ))

    fig.update_layout(
        title=f'Top {top_n} Orthogroups Counts by Species',
        xaxis=dict(title='Ortho_ID'),
        yaxis=dict(title='Counts'),
        barmode='stack',
        font=dict(size=12)
    )

    if italic_labels:
        for trace in fig.data:
            trace.name = f"<i>{trace.name}</i>"

    if config.save_plots:
        if custom_filename:
            fig.write_image(
                f"{config.output_path}/{custom_filename}.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/{custom_filename}.html")
        else:
            fig.write_image(
                f"{config.output_path}/top_{top_n}_ortho_groups.{config.file_format}", scale=config.dpi/96)
            fig.write_html(
                f"{config.output_path}/top_{top_n}_ortho_groups.html")

    if config.show:
        fig.show()


def plot_orthogroup_expression(adata: AnnData, config: PlotConfig, ortho_id: str, custom_filename: str = None,
                               trait: str = "tissue") -> None:
    """
    Plots a heatmap of the expression levels of transcripts belonging to a specific orthogroup across different tissues,
    with clustering of columns only.

    Parameters:
    - adata (AnnData): AnnData object containing the transcript expression data.
    - config (PlotConfig): Configuration object for plot control.
    - ortho_id (str): The orthogroup ID to analyze.
    - custom_filename (str): Custom filename for the plot.
    - trait (str): The trait column to use for grouping the data.
    """

    if trait not in adata.obs.columns:
        trait = "Combined_Trait"
        if trait not in adata.obs.columns:
            raise ValueError(
                f"Trait column {trait} not found in the AnnData object.")

    ortho_transcripts = adata.var[adata.var['ortho_ID'] == ortho_id].index

    # Get expression data for the orthogroup transcripts
    expression_data = adata[:, ortho_transcripts].to_df()
    expression_data[trait] = adata.obs[trait]
    expression_data = expression_data.groupby(trait).mean()

    g = sns.clustermap(expression_data.T, cmap='viridis', cbar_kws={'label': 'Expression level', 'orientation': 'vertical'},
                       col_cluster=True, row_cluster=True)
    g.figure.suptitle(
        f'Expression of Orthogroup {ortho_id} Transcripts Across Tissues', y=1.02)
    g.ax_heatmap.set_xlabel('Tissue')
    g.ax_heatmap.set_ylabel('Transcripts')

    # Adjust colorbar position to be smaller and to the right corner
    # [left, bottom, width, height]
    g.cax.set_position([0.88, 0.88, 0.02, 0.1])

    if config.save_plots:
        if custom_filename:
            g.savefig(
                f"{config.output_path}/{custom_filename}.{config.file_format}", dpi=config.dpi)
        else:
            g.savefig(
                f"{config.output_path}/expression_{ortho_id}.{config.file_format}", dpi=config.dpi)

    # Show plot if required
    if config.show:
        plt.show()

    if config.save_raw:
        raw_data_filename = custom_filename if custom_filename else f"orthogroup_expression_{ortho_id}"
        raw_data_filepath = f"{config.output_path}/{raw_data_filename}.csv"
        expression_data.to_csv(raw_data_filepath)

    plt.close()


def plot_hog_levels(adatas: List[AnnData], config: PlotConfig, custom_filename: str = None) -> None:
    """
    Plots the HOG levels per species and common HOGs using Plotly for interactive visualization.

    Parameters:
    - adatas (List[AnnData]): List of AnnData objects containing the HOG data.
    - config (PlotConfig): Configuration object for plot control.
    - custom_filename (str): Custom filename for the plot.
    """

    # Prepare data for plotting
    species_hog_levels = extract_hog_levels_per_species(adatas)
    common_hogs_per_level = get_common_hogs_per_level(species_hog_levels)

    fig = go.Figure()

    levels = list(next(iter(species_hog_levels.values())).keys())

    for species, hog_levels in species_hog_levels.items():
        num_hogs_per_level = [len(hog_levels[level]) for level in levels]
        fig.add_trace(go.Scatter(
            x=levels,
            y=num_hogs_per_level,
            mode='lines+markers',
            name=species
        ))

    common_num_hogs_per_level = [
        len(common_hogs_per_level[level]) for level in levels]
    fig.add_trace(go.Scatter(
        x=levels,
        y=common_num_hogs_per_level,
        mode='lines+markers',
        name='Common',
        line=dict(dash='dash')
    ))

    fig.update_layout(
        title='Hierarchical Orthologous Groups (HOG) Levels per Species',
        xaxis_title='Levels',
        yaxis_title='Number of HOGs',
        template='plotly_white'
    )

    if config.show:
        fig.show()

    if config.save_plots:
        if custom_filename:
            fig.write_image(
                f"{config.output_path}/{custom_filename}.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/{custom_filename}.html")
        else:
            fig.write_image(
                f"{config.output_path}/hog_levels.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/hog_levels.html")


def plot_hog_eigengenes(adatas: List[AnnData], config: PlotConfig, module: str, plant1: str, plant2: str, term1: str, term2: str, custom_filename: str = None, common: bool = True) -> None:
    """
    Plots line graph of HOG eigengenes for selected modules of two plants.

    Parameters:
    - adatas (List[AnnData]): List of AnnData objects, one per plant.
    - config (PlotConfig): Configuration object for plot control.
    - module (str): Column name in `adata.var` to group by, 'moduleColors'.
    - plant1 (str): Name of the first plant species.
    - plant2 (str): Name of the second plant species.
    - term1 (str): Module for the first plant.
    - term2 (str): Module for the second plant.
    - custom_filename (str): Custom filename for the plot.
    - common (bool): Whether to include only common HOGs between the two plants.
    """

    aggregated_df1 = None
    aggregated_df2 = None

    for adata in adatas:
        species = adata.uns['species']
        if species == plant1:
            aggregated_df1 = aggregate_eigengenes_filter(adata, module, term1)
        elif species == plant2:
            aggregated_df2 = aggregate_eigengenes_filter(adata, module, term2)

    if aggregated_df1 is None or aggregated_df2 is None:
        print(f"Error: Aggregated data for one or both plants not found.")
        return

    if common:
        common_hogs = get_common_hogs(aggregated_df1, aggregated_df2)
        aggregated_df1 = aggregated_df1.loc[common_hogs]
        aggregated_df2 = aggregated_df2.loc[common_hogs]

    df1 = pd.DataFrame({
        'HOGs': aggregated_df1.index,
        'Module Membership': aggregated_df1.mean(axis=1),
        'Plant': f"{plant1} - {term1}"
    })
    df2 = pd.DataFrame({
        'HOGs': aggregated_df2.index,
        'Module Membership': aggregated_df2.mean(axis=1),
        'Plant': f"{plant2} - {term2}"
    })

    combined_df = pd.concat([df1, df2])

    fig = px.line(combined_df, x='HOGs', y='Module Membership', color='Plant',
                  title=f'Module Membership for Module Eigengenes of {term1} and {term2} in <i>{plant1}</i> and <i>{plant2}</i>')

    if config.show:
        fig.show()

    if config.save_plots:
        if custom_filename:
            fig.write_image(
                f"{config.output_path}/{custom_filename}.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/{custom_filename}.html")
        else:
            fig.write_image(
                f"{config.output_path}/hog_eigengenes.{config.file_format}", scale=config.dpi/96)
            fig.write_html(f"{config.output_path}/hog_eigengenes.html")


def plot_cluster_expression(adata: AnnData, config: PlotConfig, cluster_id: int, cluster_name: str,
                            custom_filename: str = None, transform: str = None, trait: str = "tissue") -> None:
    """
    Plots a heatmap of the expression levels of transcripts belonging to a specific cluster across different tissues,
    with clustering of columns only.

    Parameters:
    - adata (AnnData): AnnData object containing the transcript expression data.
    - config (PlotConfig): Configuration object for plot control.
    - cluster_id (int): The cluster ID to analyze.
    - cluster_name (str): The name of the column in adata.var containing the cluster assignments.
    - custom_filename (str): Custom filename for the plot.
    - transform (str): Transformation to apply to the data ('log' or 'z').
    - trait (str): Column name in adata.obs containing the trait information.
    """

    if trait not in adata.obs.columns:
        trait = "Combined_Trait"
        if trait not in adata.obs.columns:
            raise ValueError(
                f"Trait column {trait} not found in the AnnData object.")

    cluster_transcripts = adata.var[adata.var[cluster_name]
                                    == cluster_id].index

    # Get expression data for the cluster transcripts
    expression_data = adata[:, cluster_transcripts].to_df()
    expression_data[trait] = adata.obs[trait]
    expression_data = expression_data.groupby(trait).mean()

    # Apply transformation if specified
    if transform == 'log':
        expression_data = np.log2(expression_data + 1)
        cbar_label = 'Log2(Expression level + 1)'
    elif transform == 'z':
        expression_data = zscore(expression_data, axis=1)
        cbar_label = 'Z-Score of Expression level'
    else:
        cbar_label = 'Expression level'

    g = sns.clustermap(expression_data.T, cmap='viridis', cbar_kws={'label': cbar_label, 'orientation': 'vertical'},
                       col_cluster=True, row_cluster=True)
    g.figure.suptitle(
        f'Expression of Cluster {cluster_id} Transcripts Across Tissues', y=1.02)
    g.ax_heatmap.set_xlabel('Tissue')
    g.ax_heatmap.set_ylabel('Transcripts')

    # Adjust colorbar position to be smaller and to the right corner
    # [left, bottom, width, height]
    g.cax.set_position([0.88, 0.88, 0.02, 0.1])

    if config.save_plots:
        if custom_filename:
            g.savefig(
                f"{config.output_path}/{custom_filename}.{config.file_format}", dpi=config.dpi)
        else:
            g.savefig(
                f"{config.output_path}/{cluster_name}_{cluster_id}_heatmap.{config.file_format}", dpi=config.dpi)

    # Show plot if required
    if config.show:
        plt.show()

    if config.save_raw:
        raw_data_filename = custom_filename if custom_filename else f"cluster_expression_{cluster_id}"
        raw_data_filepath = f"{config.output_path}/{raw_data_filename}.csv"
        expression_data.to_csv(raw_data_filepath)

    plt.close()

# Based on barplotModuleEigenGene function from PyWGCNA "https://github.com/mortazavilab/PyWGCNA/blob/main/PyWGCNA/wgcna.py"


def plot_eigengenes(adata: AnnData, eigengenes: pd.DataFrame, config: PlotConfig, column: str = "clusters",
                    trait: str = "tissue", custom_filename: str = None, custom_stages: List[str] = None) -> None:
    """
    Bar plot of module eigengene figure in given eigengenes DataFrame.

    Parameters:
    - adata (AnnData): AnnData object used for static information and traits.
    - eigengenes (pd.DataFrame): DataFrame containing the module eigengenes.
    - config (PlotConfig): Configuration object for plot control.
    - column (str): Column name used for plotting purposes.
    - trait (str): Column name in adata.obs containing the trait information.
    - custom_filename (str): Custom filename for the plot.
    - custom_stages (List[str]): List of custom stages to use for coloring.
    - trait (str): The trait column to use for grouping the data.
    """

    plt.rcParams['savefig.bbox'] = 'tight'

    if trait not in adata.obs.columns:
        trait = "Combined_Trait"
        if trait not in adata.obs.columns:
            print(
                f"Trait {trait} not found in adata.obs. Skipping plot_eigengenes function.")
            return

    # Use the trait information from adata.obs
    sample_info = adata.obs[trait]
    metadata_colors = generate_stage_color_dict(custom_stages=custom_stages)
    file_paths = []

    if sample_info.shape[0] == 0:
        print("WARNING: There is no sample information in the given object! Skipping plot_eigengenes function.")
        return None

    # Calculate the number of valid samples (non-null) for each module
    cluster_counts = eigengenes.notnull().sum().sort_values()

    for module_name in cluster_counts.index:
        ME = pd.DataFrame(eigengenes[module_name].values, columns=[
                          'eigengeneExp'], index=eigengenes.index)

        if module_name.startswith(f"{adata.uns['name']}: "):
            module_name = module_name[len(f"{adata.uns['name']}: "):]

        df = ME.copy(deep=True)
        df['all'] = sample_info
        ybar = df[['all', 'eigengeneExp']].groupby(['all']).mean()[
            'eigengeneExp']
        ebar = df[['all', 'eigengeneExp']].groupby(['all']).std()[
            'eigengeneExp']
        label = list(ybar.index)
        dot = df[['all', 'eigengeneExp']].copy()
        ind = {val: i for i, val in enumerate(label)}
        dot.replace(ind, inplace=True)
        xdot = dot['all']
        ydot = dot['eigengeneExp']

        if metadata_colors is None:
            palette = "lightblue"
        else:
            palette = [metadata_colors.get(val, "lightblue") for val in label]

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(len(label) + 2, 4))

        title = f"{adata.uns['name']}: Module Eigengene for {module_name.capitalize()}"
        file_name = f"{custom_filename}_barplot_eigengene_{module_name}" if custom_filename else f"barplot_eigengene_{module_name}"

        if column == "eigengenes":
            if module_name.startswith("ME"):
                module_name = module_name[2:]
            title = f"Module Eigengene for {module_name}"
            file_name = f"module_barplot_eigengene_{module_name}"
        elif column == "modules":
            title = f"Module Eigengene for {module_name}"

        ind = list(range(len(label)))
        ax.tick_params(axis='y', labelsize=15)
        ax.set_title(title, size=20, fontweight="bold")
        ax.bar(ind, ybar, align='center', color=palette)
        ax.errorbar(ind, ybar, yerr=ebar, fmt="o", color="r")
        ax.scatter(xdot, ydot, c='black', alpha=0.5)
        ax.set_xticks(np.arange(len(ind)))
        ax.set_xticklabels(label, rotation=90, fontsize=15)
        ax.set_ylabel('Eigengene Expression', fontsize=15)
        ax.set_facecolor('white')

        for spine in ax.spines.values():
            spine.set_visible(False)

        fig.subplots_adjust(bottom=0.3)

        if config.show:
            plt.show()
        if config.save_plots:
            fig.savefig(
                f"{config.output_path}/{file_name}.{config.file_format}", dpi=config.dpi)
            file_paths.append(
                f"{config.output_path}/{file_name}.{config.file_format}")

        plt.close(fig)

    return file_paths

# Based on module_trait_relationships_heatmap function from PyWGCNA "https://github.com/mortazavilab/PyWGCNA/blob/main/PyWGCNA/wgcna.py"


def plot_module_trait_relationships_heatmap(adata: AnnData, eigengenes: pd.DataFrame, module_info: pd.DataFrame,
                                            config: PlotConfig, meta_data: Union[str, List[str]], topic,
                                            alternative: str = 'two-sided', custom_filename=None) -> None:
    """
    Plot a heatmap of the correlation between module eigengenes and traits.

    Parameters:
    - adata (AnnData): AnnData object, used for static information (e.g., title from adata.uns) and traits from adata.obs.
    - eigengenes (pd.DataFrame): DataFrame containing module eigengenes.
    - module_info (pd.DataFrame): DataFrame containing module information, including total transcripts.
    - config (PlotConfig): Configuration object for plot control.
    - meta_data (str or list): Column name(s) in adata.obs containing the metadata for the trait analysis.
    - topic (str): The topic of the trait analysis.
    - alternative (str): Alternative hypothesis for the correlation test.
    - custom_filename (str): Custom filename for the plot.
    """

    plt.rcParams['savefig.bbox'] = 'tight'

    if isinstance(meta_data, str):
        meta_data = [meta_data]

    traits_df = get_data_trait_df(adata, meta_data)

    trait_corr = pd.DataFrame(index=eigengenes.columns,
                              columns=traits_df.columns, dtype="float")
    trait_pval = pd.DataFrame(index=eigengenes.columns,
                              columns=traits_df.columns, dtype="float")

    for i in eigengenes.columns:
        for j in traits_df.columns:
            tmp = stats.pearsonr(
                eigengenes[i], traits_df[j], alternative=alternative)
            trait_corr.loc[i, j] = tmp[0]
            trait_pval.loc[i, j] = tmp[1]

    num_modules = trait_pval.shape[0]
    num_traits = trait_pval.shape[1]

    figsize = (max(20, int(num_modules * 0.8)), num_traits * 0.8)

    fig, ax = plt.subplots(figsize=figsize, facecolor='white')

    # Get the counts of each cluster from module_info
    cluster_counts = module_info['total_transcripts'].to_dict()
    all_clusters_count = cluster_counts.get("All Clusters", 0)

    xlabels = []
    for label in eigengenes.columns:
        if label.lower() == 'all clusters':
            xlabels.append(f"All Clusters ({all_clusters_count})")
        else:
            xlabels.append(
                f"{label.capitalize()} ({cluster_counts.get(label, 0)})")

    ylabels = traits_df.columns

    tmp_cor = trait_corr.T.round(decimals=2)
    tmp_pvalue = trait_pval.T.round(decimals=2)
    labels = (np.asarray(["{0}\n({1})".format(cor, pvalue) for cor, pvalue in zip(
        tmp_cor.values.flatten(), tmp_pvalue.values.flatten())])).reshape(trait_corr.T.shape)

    font_scale = max(0.5, 1.5 / max(num_modules, num_traits))
    annot_font_size = min(max(10, 300 / (num_modules + num_traits)), 20)
    label_font_size = min(max(10, 300 / (num_modules + num_traits)), 20)
    title_font_size = min(max(20, 500 / (num_modules + num_traits)), 30)
    legend_font_size = min(max(10, 300 / num_traits), 20)

    sns.set_theme(font_scale=font_scale)
    res = sns.heatmap(trait_corr.T, annot=labels, fmt="", cmap='viridis',
                      vmin=-1, vmax=1, ax=ax, annot_kws={'size': annot_font_size, "weight": "bold"},
                      xticklabels=xlabels, yticklabels=ylabels)

    res.set_xticklabels(res.get_xmajorticklabels(
    ), fontsize=label_font_size, fontweight="bold", rotation=90)
    res.set_yticklabels(res.get_ymajorticklabels(),
                        fontsize=label_font_size, fontweight="bold")
    plt.yticks(rotation=0)

    ax.set_title(f"{adata.uns['name']}: Module-trait Relationships Heatmap",
                 fontsize=title_font_size, fontweight="bold")
    ax.set_facecolor('white')

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=legend_font_size)

    fig.tight_layout()

    file_name = f"{custom_filename}_trait_relationships_heatmap" if custom_filename else f"{topic}_trait_relationships_heatmap"

    full_path = f"{config.output_path}/{file_name}.{config.file_format}"

    if config.show:
        plt.show()
    if config.save_plots:
        fig.savefig(full_path)

    plt.close(fig)

    return full_path


def plot_eigengene_expression(df: pd.DataFrame, config: PlotConfig, custom_filename: str = None,
                              width: int = 1200, height: int = 800) -> Figure:
    """
    Creates an interactive Plotly figure showing the eigengene expression across tissues, species, and clusters.

    Parameters:
    - df (pd.DataFrame): The DataFrame containing the eigengene expression data.
    - config (PlotConfig): Configuration object for plot control.
    - custom_filename (str): Custom filename for the plot.

    Returns:
    - Figure: The Plotly figure object representing the interactive visualization.
    """

    # Create a combined identifier for Species and Cluster
    df['Species_Cluster'] = df['Species'] + ' - ' + df['Cluster']

    # Create the interactive Plotly figure
    fig = px.box(
        df,
        x='Tissue',
        y='Expression',
        color='Species_Cluster',
        title='Eigengene Expression Across Tissues, Species, and Clusters',
        labels={
            'Expression': 'Eigengene Expression',
            'Tissue': 'Tissue',
            'Species_Cluster': 'Species and Cluster'
        },
        hover_data=['Sample', 'Species', 'Cluster']
    )

    # Hide traces by default
    for i, trace in enumerate(fig.data):
        if i > 0:
            trace.visible = "legendonly"

    # Update layout for better readability
    fig.update_layout(
        xaxis_title='Tissue',
        yaxis_title='Eigengene Expression',
        xaxis_tickangle=45,
        legend_title='Species and Cluster',
        height=height,
        width=width,
        margin=dict(l=50, r=50, t=100, b=150),
        font=dict(size=12)
    )

    # Adjust legend
    fig.update_layout(
        legend=dict(
            title='Species and Cluster',
            orientation='v',
            xanchor='left',
            x=1.1,
            yanchor='top',
            y=1
        ),
        margin=dict(l=50, r=150, t=100, b=150)
    )

    fig.update_xaxes(showgrid=True)
    fig.update_yaxes(showgrid=True)

    custom_filename = custom_filename if custom_filename else "eigengene_expression"
    full_path = f"{config.output_path}/{custom_filename}.{config.file_format}"
    html_path = f"{config.output_path}/{custom_filename}.html"

    if config.show:
        fig.show()

    if config.save_plots:
        fig.write_image(full_path, scale=config.dpi/96)
        fig.write_html(html_path)

    if config.save_raw:
        raw_data_filepath = f"{config.output_path}/{custom_filename}.csv"
        df.to_csv(raw_data_filepath)

    return os.path.abspath(html_path)


def plot_eigengene_expression_bokeh(df: pd.DataFrame, config, custom_filename: str = None, metadata: str = "RanOmics") -> str:
    """
    Creates an interactive Bokeh figure showing the eigengene expression across tissues and clusters.
    Allows users to select any combination of tissues to display, and displays data as bar plots with error bars.

    Parameters:
    - df (pd.DataFrame): The DataFrame containing the eigengene expression data.
    - config: Configuration object for plot control.
    - custom_filename (str): Custom filename for the plot.
    """

    # Ensure 'Cluster' is treated as string
    df['Cluster'] = df['Cluster'].astype(str)

    # Drop rows with missing values
    df = df.dropna(subset=['Tissue', 'Expression',
                   'Sample', 'Species', 'Cluster'])

    # Create a combined identifier for Species and Cluster
    df['Species_Cluster'] = df['Species'] + ' - ' + df['Cluster']

    # Create a tuple for x-axis: (Tissue, Cluster)
    df['x'] = list(zip(df['Tissue'], df['Cluster']))

    # Prepare the data for bars and error bars
    grouped = df.groupby(['Tissue', 'Cluster'])
    mean_expr = grouped['Expression'].mean().reset_index()
    std_expr = grouped['Expression'].std().reset_index()

    mean_expr['x'] = list(zip(mean_expr['Tissue'], mean_expr['Cluster']))
    mean_expr['std'] = std_expr['Expression']

    # Calculate upper and lower values for error bars
    mean_expr['upper'] = mean_expr['Expression'] + mean_expr['std']
    mean_expr['lower'] = mean_expr['Expression'] - mean_expr['std']

    # Get the list of tissues and set default selection
    tissues = sorted(df['Tissue'].unique())
    default_selected_tissues = [tissues[0]]

    # Create full data sources (without filtering)
    full_bars_source = ColumnDataSource(data=mean_expr)
    full_points_source = ColumnDataSource(data=df)

    # Filter data sources based on the default selected tissues
    initial_bars_data = mean_expr[mean_expr['Tissue'].isin(
        default_selected_tissues)]
    initial_points_data = df[df['Tissue'].isin(default_selected_tissues)]

    bars_source = ColumnDataSource(data=initial_bars_data)
    points_source = ColumnDataSource(data=initial_points_data)

    # Get the unique x-factors for the selected tissues
    x_factors = sorted(initial_bars_data['x'].unique())

    # Create the figure
    p = figure(
        x_range=FactorRange(*x_factors),
        title='Eigengene Expression Across Traits and Clusters',
        toolbar_location="above",
        sizing_mode="stretch_both"
    )
    p.xaxis.major_label_orientation = 0.785  # Rotate x-axis labels by 45 degrees
    p.yaxis.axis_label = 'Eigengene Expression'

    # Add HoverTool for bars
    hover_bars = HoverTool(tooltips=[
        ('Trait', '@Tissue'),
        ('Cluster', '@Cluster'),
        ('Mean Expression', '@Expression'),
        ('Std Dev', '@std')
    ], renderers=[])
    p.add_tools(hover_bars)

    # Add HoverTool for individual points
    hover_points = HoverTool(tooltips=[
        ('Trait', '@Tissue'),
        ('Cluster', '@Cluster'),
        ('Expression', '@Expression'),
        ('Sample', '@Sample'),
        ('Species', '@Species')
    ], renderers=[])
    p.add_tools(hover_points)

    # Create a color mapping for clusters
    from bokeh.palettes import Category20, viridis

    clusters = sorted(df['Cluster'].unique())
    num_clusters = len(clusters)

    if num_clusters <= 20:
        colors = Category20[20][:num_clusters]
    else:
        colors = viridis(num_clusters)

    color_map = {cluster: colors[i % len(colors)]
                 for i, cluster in enumerate(clusters)}

    # Plot the bars
    bars = p.vbar(
        x='x',
        top='Expression',
        width=0.8,
        source=bars_source,
        fill_color=factor_cmap('Cluster', palette=colors,
                               factors=clusters),  # Korrektur hier
        line_color='black'
        # legend_field='Cluster'  # Removed to eliminate the legend
    )

    # Add error bars
    whisker = Whisker(base='x', upper='upper', lower='lower',
                      source=bars_source, line_width=2)
    whisker.upper_head.size = 10
    whisker.lower_head.size = 10
    p.add_layout(whisker)

    # Plot individual data points
    points = p.circle(
        x='x',
        y='Expression',
        source=points_source,
        size=6,
        color='black',
        alpha=0.6
    )

    # Update renderers for hover tools
    hover_bars.renderers.append(bars)
    hover_points.renderers.append(points)

    # MultiSelect for tissues
    tissue_multiselect = MultiSelect(
        title="Select Traits:", value=default_selected_tissues, options=tissues, size=8)

    # CustomJS callback for MultiSelect
    callback = CustomJS(args=dict(
        bars_source=bars_source,
        points_source=points_source,
        full_bars_source=full_bars_source,
        full_points_source=full_points_source,
        multiselect=tissue_multiselect,
        p=p
    ), code="""
        var active_tissues = multiselect.value;

        // Get the full data
        var full_bars_data = full_bars_source.data;
        var full_points_data = full_points_source.data;

        // Initialize new data objects
        var new_bars_data = {};
        var new_points_data = {};
        var new_x_factors = [];

        // Initialize new bars data
        for (var key in full_bars_data) {
            new_bars_data[key] = [];
        }

        // Filter bars data
        for (var i = 0; i < full_bars_data['Tissue'].length; i++) {
            if (active_tissues.includes(full_bars_data['Tissue'][i])) {
                for (var key in full_bars_data) {
                    new_bars_data[key].push(full_bars_data[key][i]);
                }
                var x_val = full_bars_data['x'][i];
                if (!new_x_factors.includes(x_val)) {
                    new_x_factors.push(x_val);
                }
            }
        }

        // Initialize new points data
        for (var key in full_points_data) {
            new_points_data[key] = [];
        }

        // Filter points data
        for (var i = 0; i < full_points_data['Tissue'].length; i++) {
            if (active_tissues.includes(full_points_data['Tissue'][i])) {
                for (var key in full_points_data) {
                    new_points_data[key].push(full_points_data[key][i]);
                }
            }
        }

        // Update data sources
        bars_source.data = new_bars_data;
        points_source.data = new_points_data;

        // Update x-axis factors
        p.x_range.factors = new_x_factors;

        // Emit changes
        bars_source.change.emit();
        points_source.change.emit();
        p.change.emit();
    """)
    tissue_multiselect.js_on_change('value', callback)

    # Layout and show with responsive sizing mode
    layout = column(tissue_multiselect, p, sizing_mode="stretch_both")

    custom_filename = custom_filename if custom_filename else "eigengene_expression_bokeh"

    # Save and show the plot
    if config.save_plots:
        output_file(f"{config.output_path}/{custom_filename}.html")
        save(layout)

    if config.show:
        show(layout)

    if config.save_raw:
        raw_data_filepath = f"{config.output_path}/{custom_filename}.csv"
        df.to_csv(raw_data_filepath)

    return os.path.abspath(f"{config.output_path}/{custom_filename}.html")


def plot_eigengene_modules_bokeh(adata, config, grouping_vars, custom_filename: str = None, width: int = None,
                                 height: int = None) -> str:
    """
    Creates an interactive Bokeh figure showing the eigengene expression across specified grouping variables.
    Allows users to select any grouping variable to display, and displays data as bar plots with error bars.
    Modules (Eigengenes) can be individually toggled via the legend.

    Parameters:
    - adata: AnnData object containing the eigengene expression data and metadata.
    - config: Configuration object for plot control.
    - grouping_vars: List of strings specifying the metadata columns that can be used for grouping.
    - custom_filename (str): Custom filename for the plot.
    - width (int): Width of the plot.
    - height (int): Height of the plot.

    Returns:
    - The absolute path to the saved HTML file.
    """

    # Extract eigengene expressions
    eigengene_df = adata.obsm['eigengenes'].copy()
    eigengene_df.index = adata.obs_names  # Ensure indices match

    # Extract metadata
    metadata_df = adata.obs.copy()
    metadata_df.index = adata.obs_names

    # Merge eigengene expressions with metadata
    df = eigengene_df.join(metadata_df)
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'sample'}, inplace=True)  # Correct renaming

    # Melt the DataFrame to long format
    id_vars = ['sample'] + list(metadata_df.columns)
    df = df.melt(id_vars=id_vars, var_name='Eigengene',
                 value_name='Expression')

    # Ensure that metadata columns are strings
    for col in metadata_df.columns:
        df[col] = df[col].astype(str)

    # List of possible grouping variables
    possible_groupings = grouping_vars

    # Create a Select widget for grouping variable
    grouping_select = Select(
        title="Select grouping variable:", value=grouping_vars[0], options=grouping_vars)

    # Function to prepare data based on selected grouping variable
    def prepare_data(selected_grouping):
        # Grouping variable
        grouping_var = selected_grouping

        # Drop rows with missing values in the grouping variable
        df_filtered = df.dropna(subset=[grouping_var, 'Expression'])

        # Prepare the data for bars and error bars
        grouped = df_filtered.groupby([grouping_var, 'Eigengene'])
        mean_expr = grouped['Expression'].mean().reset_index()
        std_expr = grouped['Expression'].std().reset_index()

        mean_expr['std'] = std_expr['Expression']
        mean_expr['upper'] = mean_expr['Expression'] + mean_expr['std']
        mean_expr['lower'] = mean_expr['Expression'] - mean_expr['std']

        # Prepare x-axis factors
        x_factors = sorted(df_filtered[grouping_var].unique())

        return df_filtered, mean_expr, x_factors

    # Initial data preparation
    selected_grouping = grouping_select.value
    df_filtered, mean_expr, x_factors = prepare_data(selected_grouping)

    # Create data sources
    bars_sources = {}
    whisker_sources = {}
    eigengenes = sorted(df['Eigengene'].unique())
    num_eigengenes = len(eigengenes)

    if num_eigengenes <= 20:
        colors = Category20[20][:num_eigengenes]
    else:
        colors = viridis(num_eigengenes)

    color_map = dict(zip(eigengenes, colors))

    # Create a figure
    p = figure(
        x_range=x_factors,
        title='Eigengene Expression',
        toolbar_location="above",
        sizing_mode="stretch_both",
        width=width,
        height=height,
    )
    p.xaxis.major_label_orientation = 0.785  # Rotate x-axis labels by 45 degrees
    p.yaxis.axis_label = 'Eigengene Expression'

    # Plot bars for each eigengene separately using dodge
    bars_renderers = []
    whisker_renderers = []
    # Total width is 0.8, divided by number of eigengenes
    bar_width = 0.8 / num_eigengenes

    for i, eigengene in enumerate(eigengenes):
        source = ColumnDataSource(
            mean_expr[mean_expr['Eigengene'] == eigengene])
        bars_sources[eigengene] = source

        # Calculate dodge offset
        # Center the groups by offsetting by half the total width
        dodge_offset = (-0.4 + (i + 0.5) * bar_width)

        bar = p.vbar(
            x=dodge(selected_grouping, dodge_offset, range=p.x_range),
            top='Expression',
            width=bar_width,
            source=source,
            fill_color=color_map[eigengene],
            line_color='black',
            legend_label=eigengene,
            muted_alpha=0.1,
            muted=False,
            name=eigengene,
        )
        bars_renderers.append(bar)

        # Add error bars
        whisker = Whisker(base=dodge(selected_grouping, dodge_offset, range=p.x_range),
                          upper='upper', lower='lower', source=source, line_width=2)
        whisker.upper_head.size = 10
        whisker.lower_head.size = 10
        p.add_layout(whisker)
        whisker_renderers.append(whisker)

    # Plot individual data points
    points_source = ColumnDataSource(data=df_filtered)
    points = p.circle(
        x=selected_grouping,
        y='Expression',
        source=points_source,
        size=6,
        color='black',
        alpha=0.6,
    )

    # Add HoverTool for bars
    hover_bars = HoverTool(
        tooltips=[
            (selected_grouping, f'@{selected_grouping}'),
            ('Eigengene', '@Eigengene'),
            ('Mean Expression', '@Expression{0.2f}'),
            ('Std Dev', '@std{0.2f}'),
        ],
        renderers=bars_renderers,
    )
    p.add_tools(hover_bars)

    # Add HoverTool for individual points
    hover_points = HoverTool(
        tooltips=[
            ('Sample', '@sample'),
            (selected_grouping, f'@{selected_grouping}'),
            ('Eigengene', '@Eigengene'),
            ('Expression', '@Expression{0.2f}'),
        ],
        renderers=[points],
    )
    p.add_tools(hover_points)

    # Adjust legend
    p.legend.title = 'Eigengene'
    p.legend.location = 'top_right'
    p.legend.click_policy = 'mute'

    # Callback to update data when grouping variable changes
    callback_code = """
        const grouping_var = cb_obj.value;
        const df = df_source.data;
        const eigengenes = eigengene_list;

        // Prepare new x_range factors
        const x_values = Array.from(new Set(df[grouping_var])).sort();
        p.x_range.factors = x_values;

        // Update bars and whiskers
        for (let i = 0; i < eigengenes.length; i++) {
            const eigengene = eigengenes[i];
            const source = bars_sources[eigengene];
            const new_data = { };

            // Filter and calculate mean and std dev
            const temp = {};
            for (let j = 0; j < df['Expression'].length; j++) {
                if (df['Eigengene'][j] === eigengene) {
                    const group = df[grouping_var][j];
                    if (!(group in temp)) {
                        temp[group] = [];
                    }
                    temp[group].push(df['Expression'][j]);
                }
            }

            const groups = Object.keys(temp);
            const means = [];
            const stds = [];
            const uppers = [];
            const lowers = [];
            const group_names = [];

            groups.forEach(group => {
                const values = temp[group];
                const mean = values.reduce((a, b) => a + b, 0) / values.length;
                const variance = values.map(x => Math.pow(x - mean, 2)).reduce((a, b) => a + b, 0) / values.length;
                const std = Math.sqrt(variance);
                means.push(mean);
                stds.push(std);
                uppers.push(mean + std);
                lowers.push(mean - std);
                group_names.push(group);
            });

            new_data[grouping_var] = group_names;
            new_data['Expression'] = means;
            new_data['std'] = stds;
            new_data['upper'] = uppers;
            new_data['lower'] = lowers;
            new_data['Eigengene'] = Array(group_names.length).fill(eigengene);

            // Update data source
            source.data = new_data;
            source.change.emit();
        }

        // Update points
        const points_data = points_source.data;
        points_data[grouping_var] = df[grouping_var];
        points_source.change.emit();

        // Update hover tooltips
        hover_bars.tooltips = [
            [grouping_var, `@${grouping_var}`],
            ['Eigengene', '@Eigengene'],
            ['Mean Expression', '@Expression{0.2f}'],
            ['Std Dev', '@std{0.2f}'],
        ];
        hover_points.tooltips = [
            ['Sample', '@sample'],
            [grouping_var, `@${grouping_var}`],
            ['Eigengene', '@Eigengene'],
            ['Expression', '@Expression{0.2f}'],
        ];

        // Update axis label
        p.xaxis.axis_label = grouping_var;
    """

    # Prepare data sources for JavaScript callback
    df_source = ColumnDataSource(df)
    bars_sources_js = {
        eigengene: bars_sources[eigengene] for eigengene in eigengenes}

    callback = CustomJS(
        args={
            'grouping_select': grouping_select,
            'p': p,
            'df_source': df_source,
            'points_source': points_source,
            'bars_sources': bars_sources_js,
            'eigengene_list': eigengenes,
            'hover_bars': hover_bars,
            'hover_points': hover_points,
        },
        code=callback_code,
    )

    grouping_select.js_on_change('value', callback)

    # Layout and display with responsive sizing mode
    layout = column(grouping_select, p, sizing_mode="stretch_both")

    custom_filename = custom_filename if custom_filename else "eigengene_expression_bokeh"

    # Save and show the plot
    if config.save_plots:
        output_file(f"{config.output_path}/{custom_filename}.html")
        save(layout)

    if config.show:
        show(layout)

    if config.save_raw:
        raw_data_filepath = f"{config.output_path}/{custom_filename}.csv"
        df.to_csv(raw_data_filepath, index=False)

    return os.path.abspath(f"{config.output_path}/{custom_filename}.html")


"""
Cytoscape.js Network Plotting
"""

NodeData = Dict[str, Union[str, int]]
Node = Dict[str, NodeData]

EdgeData = Dict[str, Union[str, float]]
Edge = Dict[str, EdgeData]


def plot_cyto_network(config: PlotConfig, custom_filename: str = "cyto_network", network_data: Dict[str, List[Union[Node, Edge]]] = None,
                      searchpath: str = "../flask/app_dev/templates/", use_colors: bool = True, use_shapes: bool = False,
                      cluster_info: Dict[str, str] = None, node_size: int = 10, highlight: List[str] = None,
                      edge_width: int = 1, highlight_color: str = "magenta",
                      use_edge_transparency: bool = False, use_edge_coloring: bool = True, filter_edges: bool = True) -> str:
    """
    Plot a network using Cytoscape.js and save it as an HTML file.

    Parameters:
    - config (PlotConfig): Configuration object for plot control.
    - custom_filename (str): Custom filename for the plot.
    - network_data (Dict[str, List[Union[Node, Edge]]]): The network data to be plotted.
    - searchpath (str): The search path for the Jinja2 template.
    - use_colors (bool): Whether to use module colors in the network plot.
    - use_shapes (bool): Whether to use different shapes for different node types.
    - cluster_info (Dict[str, str]): Dictionary containing cluster information.
    - node_size (int): The size of the nodes in the network plot.
    - edge_width (int): The width of the edges in the network plot.
    - highlight (List[str]): List of nodes (e.g. transcripts) to highlight in the network.
    - highlight_color (str): Color to use for highlighting the nodes.
    - use_edge_transparency (bool): Whether to apply transparency to edges based on weights.
    - use_edge_coloring (bool): Whether to color edges based on weights.
    - filter_edges (bool): Whether to filter edges.

    Returns:
    - The absolute path to the saved HTML file.
    """

    # If highlight is a string, convert it to a list
    if isinstance(highlight, str):
        highlight = [highlight]

    # Adding cluster information to each node
    if cluster_info:
        for node in network_data['nodes']:
            node_id = node['data']['id']
            node['data']['cluster'] = cluster_info[node_id]

    if use_shapes:
        shape_dict = generate_cyto_shape_dict(
            [node['data']['moduleColor'] for node in network_data['nodes']], None)
        for node in network_data['nodes']:
            module_color = node['data']['moduleColor']
            node['data']['shape'] = shape_dict.get(module_color, 'ellipse')

    if highlight:
        for node in network_data['nodes']:
            if node['data']['gene'] in highlight:
                node['data']['highlighted'] = True
            else:
                node['data']['highlighted'] = False

    # Load the Jinja2 environment
    env = Environment(loader=FileSystemLoader(searchpath=searchpath))

    # Render the CSS template
    css_template = env.get_template("network.css")
    css_content = css_template.render()

    # Render the JS template
    js_template = env.get_template("network.js")
    js_content = js_template.render(network_data=network_data,
                                    filter_edges=filter_edges,
                                    use_background_color=use_colors,
                                    use_shapes=use_shapes,
                                    use_cluster_tooltip=cluster_info is not None,
                                    node_size=node_size,
                                    edge_width=edge_width,
                                    highlight_color=highlight_color,
                                    use_edge_transparency=use_edge_transparency,
                                    use_edge_coloring=use_edge_coloring)

    # Render the HTML template and embed CSS and JS contents
    html_template = env.get_template("network.html")
    rendered_html = html_template.render(css_content=css_content,
                                         js_content=js_content)

    html_file_path = f"{config.output_path}/{custom_filename}.html"
    if config.show:
        display(IFrame(src=html_file_path, width='100%', height='600px'))

    if config.save_plots:
        with open((html_file_path), 'w') as f:
            f.write(rendered_html)

    if config.save_raw:
        raw_data_filename = custom_filename if custom_filename else "cyto_network"
        raw_data_filepath = f"{config.output_path}/{raw_data_filename}.json"
        with open(raw_data_filepath, 'w') as f:
            json.dump(network_data, f)

    return os.path.abspath(html_file_path)
