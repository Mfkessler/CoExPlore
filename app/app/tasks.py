# tasks.py

import os
import logging
import shutil
import wgcna.plotting as rplot
import wgcna.utils as rutils
import wgcna.ortho as rortho
import wgcna.wrapper as rwrap
from datetime import datetime
from config import Config
from flask import current_app
from .factory import celery
from .cache import get_adata_cache
from .extensions import redis_client

# Configure logger
logger = logging.getLogger(__name__)

# Plotting configuration
plot_config = rplot.PlotConfig(save_plots=True, 
                          output_path=Config.OUTPUT_DIR, 
                          file_format='png', 
                          show=False,
                          save_raw=True, 
                          dpi=300)

@celery.task(bind=True, name='init_adata_cache_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def init_adata_cache_task(self):
    logging.info("Starting init_adata_cache_task...")
    adata_cache = get_adata_cache()
    
    if adata_cache:
        logging.info("AnnData cache successfully initialized.")
        return {"status": "SUCCESS", "result": {"status": "success", "message": "AnnData cache initialized"}}
    else:
        logging.error("Failed to initialize AnnData cache.")
        return {"status": "FAILURE", "result": {"status": "error", "message": "Failed to initialize AnnData cache"}}

######################
"""Browser Analysis"""
######################

@celery.task(bind=True, name='plot_co_expression_network_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_co_expression_network_task(self, data):
    adata_cache = get_adata_cache()

    def progress_callback(message, **kwargs):
        """Helper function to send progress updates"""
        self.update_state(state='PROGRESS', meta={'status': message, **kwargs})

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    plant = data['plant']
    tom_path = f"{Config.DATA_DIR}/tom/tom_matrix_{plant}.h5"

    transcripts = data.get('transcripts')
    total_transcripts = sum(len(t_list) for t_list in transcripts.values())  # Sum of all transcript lists
    logger.info(f"Total transcripts before filtering: {total_transcripts}")

    threshold = float(data.get('threshold', 0.1))
    use_shapes = data.get('useShapesSpecies', False)
    use_colors = data.get('useColors', False)
    highlight_list = data.get('highlightList', None)
    custom_filename = f"{prefix}co_expression_network_{plant}_{threshold}"
    template_path = os.path.join(current_app.root_path, 'templates')
    topic = ""
    max_neighbors = data.get('maxNeighbors', 0)
    force_detailed_view = data.get('forceDetailedView', False)
    min_cluster_size = data.get('minClusterSize', 10)

    if force_detailed_view:
        detail_only_nodes = None
    else:
        detail_only_nodes = 500

    if max_neighbors > 0:
        logger.info(f"Max neighbors: {max_neighbors}")
        include_neighbors = True
    else:
        include_neighbors = False

    if isinstance(plant, list):
        custom_filename = f"{prefix}co_expression_network_{'_'.join(plant)}_{threshold}"

    if plant in adata_cache.keys() or isinstance(plant, list):
        if isinstance(plant, list):
            adata = [adata_cache.get_adata(p) for p in plant]
        else:
            adata = adata_cache.get_adata(plant)

        total_transcripts = sum(len(v) for v in transcripts.values())
        if total_transcripts > 8000:
            return {"status": "FAILURE", "result": {"status": "error", "message": "Too many transcripts selected, please select fewer than 8000"}}

        html_path = rwrap.analyze_co_expression_network(adata, plot_config, transcripts=transcripts, threshold=threshold, node_threshold=min_cluster_size,
                                                        obo_path=f"{Config.DATA_DIR}/go-basic.obo", topic=topic, plot_go_enrichment=False,
                                                        template_path=template_path, highlight=highlight_list, custom_filename=custom_filename,
                                                        use_colors=use_colors, use_shapes=use_shapes, progress_callback=progress_callback,
                                                        tom_prefix=f"{Config.DATA_DIR}/tom", filter_edges=False, max_neighbors=max_neighbors,
                                                        include_neighbors=include_neighbors, detail_only_nodes=detail_only_nodes)
        
        if isinstance(html_path, dict):
            return {"status": "FAILURE", "result": {"status": "error", "message": f"{html_path.get('message')}"}}

        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{plot_config.output_path}/{html_path}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded", "info": f"{adata_cache.keys()}"}}    


@celery.task(bind=True, name='plot_go_terms_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_go_terms_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    plant = data['plant']
    transcripts = data.get('transcripts')
    max_pval = float(data.get('maxPval', 0.5))
    min_fe = float(data.get('minFe', 0.2))
    min_depth = int(data.get('minDepth', 2))
    top_percentage = float(data.get('nTopPercent', 1.0))

    if plant in adata_cache.keys():
        adata = adata_cache.get_adata(plant)
        custom_filename=f"{prefix}{plant}_{max_pval}_{min_fe}_{min_depth}"
        res = rplot.plot_goea(adata, plot_config, transcripts, name=plant, custom_filename=f"{custom_filename}_tree",
                              top_percentage=top_percentage, obo_path=f"{Config.DATA_DIR}/go-basic.obo")
        
        rplot.plot_filtered_go_terms(res, plot_config, column=None, min_depth=min_depth, 
                                                 min_fold_enrichment=min_fe, p_value_threshold=max_pval,
                                                 custom_filename=custom_filename)
        output_file = f"{plot_config.output_path}/{custom_filename}_GO_Terms.html"
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}
    

@celery.task(bind=True, name='plot_venn_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_venn_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    plants = data.get('plant', [])
    transcripts = data.get('transcripts')
    adata_list = [adata_cache.get_adata(plant) for plant in plants]

    ortho_df = rortho.get_ortho_df(adata_list, transcripts)

    if adata_list and ortho_df is not None:
        logger.info(f"Transcripts: {len(transcripts)}")
        logger.info(f"Ortho_df: {ortho_df.shape}")

        species_list = [col for col in ortho_df.columns if col != 'Count']

        try:
            rplot.plot_venn(ortho_df, plot_config, species_list, custom_filename=f"{prefix}venn_diagram_{'_'.join(plants)}",
                            colors=['red', 'green', 'blue'], alpha=0.7)
            plot_filename = f"{prefix}venn_diagram_{'_'.join(plants)}.png"
            return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{plot_config.output_path}/{plot_filename}"}}
        except Exception as e:
            logger.error(f"Error in plot_venn_task: {e}")
            return {"status": "FAILURE", "result": {"status": "error", "message": "Failed to generate Venn diagram"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}
    

@celery.task(bind=True, name='plot_upset_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_upset_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    plants = data.get('plant', [])
    transcripts = data.get('transcripts')
    adata_list = [adata_cache.get_adata(plant) for plant in plants]

    ortho_df = rortho.get_ortho_df(adata_list, transcripts=transcripts)

    if adata_list and ortho_df is not None:
        logger.info(f"Transcripts: {len(transcripts)}")
        logger.info(f"Ortho_df: {ortho_df.shape}")

        species_list = [col for col in ortho_df.columns if col != 'Count']

        try:
            rplot.plot_upset(ortho_df, plot_config, species_list, custom_filename=f"{prefix}upset_plot_{'_'.join(plants)}")

            plot_filename = f"{prefix}upset_plot_{'_'.join(plants)}.png"
            return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{plot_config.output_path}/{plot_filename}"}}
        except Exception as e:
            logger.error(f"Error in plot_upset_task: {e}")
            return {"status": "FAILURE", "result": {"status": "error", "message": "Failed to generate UpSet plot"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}
    

@celery.task(bind=True, name='plot_filtered_jaccard_modules')
def plot_filtered_jaccard_modules_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    threshold = float(data.get('threshold', 0.1))
    use_shapes = data.get('useShapesSpecies', False)
    use_colors = data.get('useColors', False)
    inter_species_only = data.get('interSpeciesOnly', True)
    min_orthos = int(data.get('minOrthos', 0))

    plants = data.get('plant', [])
    if isinstance(plants, str):
        plants = [plants]
        
    transcripts = data.get('transcripts')
    adata_list = [adata_cache.get_adata(plant) for plant in plants]

    if len(adata_list) > 0:
        mds = [rutils.create_mapping_dict(adata, value_column='ortho_id', transcripts=transcripts, min_count=min_orthos) for adata in adata_list]
        j_df_modules = rutils.get_jaccard_df(mds)

        custom_filename = f"{prefix}modules_network_filtered_{threshold}"
        output_file = f"{plot_config.output_path}/{custom_filename}.html"
        color_map = "Modules" if use_colors else ""

        rplot.plot_correlation_network(j_df_modules, plot_config, color_map=color_map, threshold=threshold, 
                                       inter_species_only=inter_species_only, width=None, height=None,
                                       suffix="modules", weight_label="Jaccard Similarity", use_shapes=use_shapes,
                                       show_symbol_legend=True, show_edge_legend=True, title="Similarity of Modules between Species",
                                       custom_filename=custom_filename)
        
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}


@celery.task(bind=True, name='plot_filtered_jaccard_species')
def plot_filtered_jaccard_species_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    threshold = float(data.get('threshold', 0.1))
    use_shapes = data.get('useShapesSpecies', False)
    min_orthos = int(data.get('minOrthos', 0))

    plants = data.get('plant', [])
    if isinstance(plants, str):
        plants = [plants]
        
    transcripts = data.get('transcripts')
    adata_list = [adata_cache.get_adata(plant) for plant in plants]

    if len(adata_list) > 1:
        mds = [rutils.create_mapping_dict(adata, value_column='ortho_id', transcripts=transcripts, mapping_column="species", min_count=min_orthos) for adata in adata_list]
        j_df_species = rutils.get_jaccard_df(mds)

        custom_filename = f"{prefix}species_network_filtered_{threshold}"
        output_file = f"{plot_config.output_path}/{custom_filename}.html"
        color_map = ""

        rplot.plot_correlation_network(j_df_species, plot_config, color_map=color_map, threshold=threshold, 
                                       inter_species_only=False, width=None, height=None,
                                       suffix="species", weight_label="Jaccard Similarity", use_shapes=use_shapes,
                                       show_symbol_legend=True, show_edge_legend=True, title="Similarity of Species",
                                       custom_filename=custom_filename)
        
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}

#####################
"""Single Analysis"""
#####################

# TODO top percentage textfield in frontend
@celery.task(bind=True, name='plot_module_goea_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_module_goea_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    plant = data['plant']
    max_pval = float(data.get('maxPval', 0.5))
    min_fe = float(data.get('minFe', 0.2))
    min_depth = int(data.get('minDepth', 2))  
    top_percentage = float(data.get('nTopPercent', 1.0))

    if plant in adata_cache.keys():
        adata = adata_cache.get_adata(plant)
        column = "Module"
        output_file_base = f"{prefix}{plant}_modules_goea_{top_percentage}"
        output_file = f"{plot_config.output_path}/{output_file_base}_{column}.html"

        if f"goea_modules_{top_percentage}" not in adata.uns:
            goea_modules_key = f"goea_modules_{top_percentage}"
            rutils.add_goea_to_anndata(adata, column="module_colors", results_dir=f"{plot_config.output_path}", 
                                       uns_key=goea_modules_key, top_percentage=top_percentage, obo_path=f"{Config.DATA_DIR}/go-basic.obo")
        df = rutils.get_goea_df(adata, key=f"goea_modules_{top_percentage}", column_name="Module")

        df['P-Value'] = df['P-Value'].astype(float)
        df['Fold Enrichment'] = df['Fold Enrichment'].astype(float)
        df['Depth'] = df['Depth'].astype(int)

        rplot.plot_filtered_go_terms(df, plot_config, p_value_threshold=max_pval, min_fold_enrichment=min_fe, min_depth=min_depth, 
                                     column="Module", custom_filename=output_file_base)
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}
    

@celery.task(bind=True, name='plot_correlation_network_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_correlation_network_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    plant = data['plant']
    threshold = float(data.get('threshold', 0.5))
    use_shapes = data.get('useShapes', False)
    use_colors = data.get('useColors', False)

    color_map = "Modules" if use_colors else ""

    if plant in adata_cache.keys():
        adata = adata_cache.get_adata(plant)
        correlation_matrix = adata.obsm['eigengenes'].corr()
        custom_filename = f"{prefix}correlation_network_{plant}_{threshold}"
        output_file = f"{plot_config.output_path}/{custom_filename}.html"
        rplot.plot_correlation_network(
            correlation_matrix, plot_config, threshold=threshold,
            suffix=plant, use_shapes=use_shapes, weight_label="Eigengene Correlation",
            width=None, height=None, show_symbol_legend=False, show_edge_legend=True,
            custom_filename=custom_filename, color_map=color_map
        )
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}
    

@celery.task(bind=True, name='plot_module_props_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_module_props_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    plant = data['plant']

    if plant in adata_cache.keys():
        adata = adata_cache.get_adata(plant)
        custom_filename = f"{prefix}module_gene_proportions_{plant}"
        output_file = f"{plot_config.output_path}/{custom_filename}.html"
        rplot.plot_module_gene_proportions(adata, plot_config, custom_filename=custom_filename)
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}
    

@celery.task(bind=True, name='plot_module_orthos_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_module_orthos_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    plant = data['plant']

    if plant in adata_cache.keys():
        adata = adata_cache.get_adata(plant)
        custom_filename = f"{prefix}module_orthos_{plant}"
        output_file = f"{plot_config.output_path}/{custom_filename}.html"
        rplot.plot_orthos_and_transcripts_per_module(adata, plot_config, custom_filename=custom_filename)
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}
    

@celery.task(bind=True, name='plot_module_hubs_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_module_hubs_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    plant = data['plant']

    if plant in adata_cache.keys():
        adata = adata_cache.get_adata(plant)
        custom_filename = f"{prefix}module_hubs_{plant}"
        output_file = f"{plot_config.output_path}/{custom_filename}.html"
        rplot.plot_hub_connectivity(adata, plot_config, custom_filename=custom_filename)
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}
    
##########################
"""Comparative Analysis"""
##########################
    
@celery.task(bind=True, name='plot_jaccard_tissues_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_jaccard_tissues_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    threshold = float(data.get('threshold', 0.1))
    nTop = data.get('nTop', 100)
    use_shapes = data.get('useShapes', False)
    # use_colors = data.get('useColors', False)
    inter_species_only = data.get('interSpeciesOnly', True)

    plants = data.get('plant', [])
    adata_list = [adata_cache.get_adata(plant) for plant in plants]

    if len(adata_list) > 1:
        dicts = rutils.get_top_expressed_genes_by_tissue(adata_list, n=nTop, column='ortho_id')
        j_df_tissues = rutils.get_jaccard_df(dicts)

        custom_filename = f"{prefix}tissues_network_{threshold}"
        output_file = f"{plot_config.output_path}/{custom_filename}.html"
        rplot.plot_correlation_network(j_df_tissues, plot_config, threshold=threshold, inter_species_only=inter_species_only, 
                                       weight_label="Jaccard Similarity", width=None, height=None,
                                       custom_filename=custom_filename, use_shapes=use_shapes, show_symbol_legend=True, 
                                       show_edge_legend=True, title="Similarity of Tissues between Species")
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}


@celery.task(bind=True, name='plot_tissues_corr_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_tissues_corr_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    nTop = data.get('nTop', 100)
    column = "Tissues"

    plants = data.get('plant', [])
    adata_list = [adata_cache.get_adata(plant) for plant in plants]

    if len(adata_list) > 1:
        dicts = rutils.get_top_expressed_genes_by_tissue(adata_list, n=nTop, column='ortho_id')
        j_df_tissues = rutils.get_jaccard_df(dicts)

        custom_filename = f"{prefix}tissues_heatmap_{nTop}"
        output_file = f"{plot_config.output_path}/{custom_filename}.html"
        rplot.plot_overlap(j_df_tissues, plot_config, column=column, custom_filename=custom_filename, title="Similarity of Tissues between Species",
                           x_label="Tissues", y_label="Tissues", height=None, width=None)
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}
    

@celery.task(bind=True, name='plot_jaccard_modules_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_jaccard_modules_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    threshold = float(data.get('threshold', 0.1))
    use_shapes = data.get('useShapes', False)
    use_colors = data.get('useColors', False)
    inter_species_only = data.get('interSpeciesOnly', True)

    plants = data.get('plant', [])
    adata_list = [adata_cache.get_adata(plant) for plant in plants]

    if len(adata_list) > 1:
        mds = [rutils.create_mapping_dict(adata, value_column='ortho_id') for adata in adata_list]
        j_df_modules = rutils.get_jaccard_df(mds)

        custom_filename = f"{prefix}modules_network_{threshold}"
        output_file = f"{plot_config.output_path}/{custom_filename}.html"
        color_map = "Modules" if use_colors else ""

        rplot.plot_correlation_network(j_df_modules, plot_config, color_map=color_map, threshold=threshold, 
                                       inter_species_only=inter_species_only, width=None, height=None,
                                       suffix="modules", weight_label="Jaccard Similarity", use_shapes=use_shapes,
                                       show_symbol_legend=True, show_edge_legend=True, title="Similarity of Modules between Species",
                                       custom_filename=custom_filename)
        
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}
    

@celery.task(bind=True, name='plot_modules_corr_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_modules_corr_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    column = "Modules"

    plants = data.get('plant', [])
    adata_list = [adata_cache.get_adata(plant) for plant in plants]

    if len(adata_list) > 1:
        mds = [rutils.create_mapping_dict(adata, value_column='ortho_id') for adata in adata_list]
        j_df_modules = rutils.get_jaccard_df(mds)

        custom_filename = f"{prefix}modules_heatmap"
        output_file = f"{plot_config.output_path}/{custom_filename}.html"
        rplot.plot_overlap(j_df_modules, plot_config, column=column, custom_filename=custom_filename, title="Similarity of Modules between Species",
                           x_label="Modules", y_label="Modules", height=None, width=None)
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}


@celery.task(bind=True, name='plot_hog_level_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_hog_level_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    plants = data.get('plant', [])
    adata_list = [adata_cache.get_adata(plant) for plant in plants]

    if len(adata_list) > 1:
        custom_filename = f"{prefix}hog_levels"
        output_file = f"{plot_config.output_path}/{custom_filename}.html"
        rplot.plot_hog_levels(adata_list, plot_config, custom_filename=custom_filename)
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}
    

@celery.task(bind=True, name='plot_species_correlation_network_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def plot_species_correlation_network_task(self, data):
    adata_cache = get_adata_cache()

    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    plot_config.output_path = user_dir
    prefix = data.get('prefix', '')

    threshold = float(data.get('threshold', 0.5))
    pval_threshold = float(data.get('maxPval', 0.05))
    positive = data.get('positive', True)
    use_shapes = data.get('useShapes', False)
    use_colors = data.get('useColors', False)
    inter_species_only = data.get('interSpeciesOnly', True)

    plants = data.get('plant', [])
    adata_list = [adata_cache.get_adata(plant) for plant in plants]

    comp_dir = f"{Config.DATA_DIR}/compare_eigengenes"

    if use_colors:
        color_map = "Modules"
    else:
        color_map = ""

    if len(adata_list) > 1:
        corr, pvals, combined_df = rutils.load_comparison(adata_list, comp_dir)

        custom_filename = f"{prefix}module_membership_network_{threshold}"
        output_file = f"{plot_config.output_path}/{custom_filename}.html"
        rplot.plot_correlation_network(corr, plot_config, threshold=threshold, pvals=pvals, pval_threshold=pval_threshold, 
                                       inter_species_only=inter_species_only, positive=positive, color_map=color_map, suffix="species",
                                       show_symbol_legend=True, show_edge_legend=True, custom_filename=custom_filename,
                                       weight_label="Correlation", use_shapes=use_shapes, width=None, height=None,
                                       title="Module Membership Correlation of Modules between Species")
        return {"status": "SUCCESS", "result": {"status": "success", "plot_url": f"{Config.BASE_URL}/{output_file}"}}
    else:
        return {"status": "FAILURE", "result": {"status": "error", "message": "Data not loaded"}}

########################
"""Session Management"""
########################

# Cleanup old sessions using redis
@celery.task(bind=True, name='cleanup_old_sessions_task', queue=Config.CELERY_TASK_DEFAULT_QUEUE)
def cleanup_old_sessions_task(self):
    print("Cleaning up old sessions")
    current_time = datetime.now().timestamp()  # Current time as a timestamp

    last_heartbeat = redis_client.hgetall('last_heartbeat')  # Get all session heartbeats

    for session_id, timestamp in last_heartbeat.items():
        session_id = session_id.decode('utf-8')
        timestamp = float(timestamp.decode('utf-8'))  # Convert the timestamp back to float

        session_age = current_time - timestamp
        if session_age > 60 * 60:  # If the session is older than 1 hour
            user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
            print(f"Attempting to remove directory: {user_dir}")
            
            if os.path.exists(user_dir):
                try:
                    shutil.rmtree(user_dir)
                    print(f"Directory {user_dir} removed due to inactivity")
                except Exception as e:
                    print(f"Error removing directory {user_dir}: {e}")
            else:
                print(f"Directory {user_dir} does not exist")

            # Remove the session from Redis
            redis_client.hdel('last_heartbeat', session_id)
        else:
            print(f"Session {session_id} is still active")