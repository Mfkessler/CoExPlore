import pytest
import os
import logging
from anndata import read_h5ad
from wgcna.wrapper import analyze_co_expression_network, get_tom_data
from wgcna.plotting import PlotConfig

logging.basicConfig(level=logging.INFO)

# Test data paths
BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")

script_dir = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.abspath(os.path.join(script_dir, "../../app/app/templates"))

if not os.path.exists(TEMPLATE_DIR):
    raise FileNotFoundError(f"Template directory not found: {TEMPLATE_DIR}")

@pytest.fixture
def ann_data_list():
    """Fixture: Load multiple AnnData objects."""
    # Load all AnnData objects from BASE_DIR/data (Folders S1 to S2) using a list comprehension
    adatas = [read_h5ad(os.path.join(BASE_DIR, "data", f"S{i}", f"S{i}.h5ad")) for i in range(1, 3)]

    return adatas

@pytest.fixture
def ann_data_single():
    """Fixture: Load a single AnnData object."""
    return read_h5ad(os.path.join(BASE_DIR, "data", "S1", "S1.h5ad"))

@pytest.fixture
def plot_config():
    """Fixture: Standard plot configuration."""
    output_path = os.path.join(BASE_DIR, "output")
    os.makedirs(output_path, exist_ok=True)  # Ensure the folder exists

    return PlotConfig(output_path=output_path, show=False, save_plots=True)

# Tests

def test_get_tom_data_single(ann_data_single):
    """Test get_tom_data with a single AnnData object."""
    tom = get_tom_data(tom_path=os.path.join(BASE_DIR, "data", "S1", "tom_matrix.h5"), 
                       adata=ann_data_single, 
                       query="",
                       threshold=0.2)
    
    assert tom is not None, "TOM should be loaded."
    assert tom.shape[0] == tom.shape[1], "TOM should be square."

def test_get_tom_data_list(ann_data_list):
    """Test get_tom_data with a list of AnnData objects."""
    toms, adatas = get_tom_data(tom_path=[os.path.join(BASE_DIR, "data", f"S{i}", "tom_matrix.h5") for i in range(1, 9)], 
                                adata=ann_data_list, 
                                query="",
                                threshold=0.2)
    
    assert len(toms) == len(ann_data_list), "There should be one TOM per AnnData object."
    
    for tom in toms:
        assert tom.shape[0] == tom.shape[1], "TOM should be square."

def test_analyze_co_expression_network_single(ann_data_single, plot_config):
    """Test analyze_co_expression_network with a single AnnData object."""
    logging.info(f"PlotConfig: {plot_config.output_path}")
    result = analyze_co_expression_network(
        adata=ann_data_single,
        config=plot_config,
        topic="Test1",
        tom_path=os.path.join(BASE_DIR, "data", "S1", "tom_matrix.h5"),
        out="html",
        template_path=TEMPLATE_DIR,
        plot_go_enrichment=False,
        query=""
    )

    assert isinstance(result, str), "The output should be a file path."
    assert result.endswith(".html"), "The output should be an HTML file."

def test_analyze_co_expression_network_list(ann_data_list, plot_config):
    """Test analyze_co_expression_network with a list of AnnData objects."""
    result = analyze_co_expression_network(
        adata=ann_data_list,
        config=plot_config,
        topic="Test2",
        tom_path=[os.path.join(BASE_DIR, "data", f"S{i}", "tom_matrix.h5") for i in range(1, 9)],
        out="html",
        template_path=TEMPLATE_DIR,
        plot_go_enrichment=False,
        query=""
    )

    assert isinstance(result, str), "The output should be a file path."
    assert result.endswith(".html"), "The output should be an HTML file."

def test_invalid_tool(ann_data_single, plot_config):
    """Test analyze_co_expression_network with an invalid tool."""
    with pytest.raises(ValueError, match="Invalid tool"):
        analyze_co_expression_network(
            adata=ann_data_single,
            config=plot_config,
            topic="Test Topic",
            tom_path=os.path.join(BASE_DIR, "data", "S1", "tom_matrix.h5"),
            tool="invalid_tool"
        )
