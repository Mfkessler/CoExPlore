# task_mapping.py

from .tasks import (
    plot_venn_task,
    plot_upset_task,
    plot_filtered_jaccard_modules_task,
    plot_filtered_jaccard_species_task,
    plot_module_goea_task,
    plot_correlation_network_task,
    plot_module_props_task,
    plot_module_orthos_task,
    plot_module_hubs_task,
    plot_jaccard_tissues_task,
    plot_tissues_corr_task,
    plot_jaccard_modules_task,
    plot_modules_corr_task,
    plot_hog_level_task,
    plot_species_correlation_network_task,
    plot_co_expression_network_task,
    plot_go_terms_task
)

# Map all analysis types to their respective tasks
analysis_tasks = {
    'plot_venn': plot_venn_task,
    'plot_upset': plot_upset_task,
    'plot_filtered_jaccard_modules': plot_filtered_jaccard_modules_task,
    'plot_filtered_jaccard_species': plot_filtered_jaccard_species_task,
    'plot_module_goea': plot_module_goea_task,
    'plot_correlation_network': plot_correlation_network_task,
    'plot_module_props': plot_module_props_task,
    'plot_module_orthos': plot_module_orthos_task,
    'plot_module_hubs': plot_module_hubs_task,
    'plot_jaccard_tissues': plot_jaccard_tissues_task,
    'plot_tissues_corr': plot_tissues_corr_task,
    'plot_jaccard_modules': plot_jaccard_modules_task,
    'plot_modules_corr': plot_modules_corr_task,
    'plot_hog_level': plot_hog_level_task,
    'plot_species_correlation_network': plot_species_correlation_network_task,
    'plot_co_expression_network': plot_co_expression_network_task,
    'plot_go_terms': plot_go_terms_task
}
