import pytest
import pandas as pd
from unittest.mock import mock_open, patch
import wgcna.utils as rurils

# ## Test cases for find_species_by_initials
# Sample content that mimics the species file structure
sample_content = """Homo sapiens
Pan troglodytes
Mus musculus"""


def test_find_species_valid_initials():
    # Test case for valid initials
    with patch("builtins.open", mock_open(read_data=sample_content)):
        assert rurils.find_species_by_initials("Hs") == "Homo sapiens"
        assert rurils.find_species_by_initials("Pt") == "Pan troglodytes"


def test_find_species_no_match():
    # Test case for initials that do not match
    with patch("builtins.open", mock_open(read_data=sample_content)):
        with pytest.raises(ValueError) as excinfo:
            rurils.find_species_by_initials("Xx")
            
        assert "No species found with the provided initials" in str(excinfo.value)


def test_find_species_wrong_format():
    # Test case for input format errors
    result = rurils.find_species_by_initials("Homo")

    assert result == "Homo", "Function should return the input if format is invalid"


def test_find_species_empty_file():
    # Test case for an empty file
    with patch("builtins.open", mock_open(read_data="")):
        with pytest.raises(ValueError) as excinfo:
            rurils.find_species_by_initials("Hs")

        assert "No species found with the provided initials" in str(excinfo.value)


def test_file_not_found():
    # Test case when the file is not found
    with patch("builtins.open", side_effect=FileNotFoundError):
        with pytest.raises(FileNotFoundError):
            rurils.find_species_by_initials("Hs")


# ## Test cases for map_tissue_types
def test_correct_mapping():
    # Create a DataFrame with sample data
    df = pd.DataFrame({"sample": ["X_A_Y", "X_B_Y", "X_C_Y"]})
    # Run the function
    result = rurils.map_tissue_types(df)

    assert list(result["tissue"]) == ["Bud stage 1", "Bud stage 2", "Bud stage 3"]
    assert "sample" in result.columns

def test_no_match():
    df = pd.DataFrame({"sample": ["X_Z_Y"]})
    result = rurils.map_tissue_types(df)

    assert list(result["tissue"]) == ["Unknown"]

def test_empty_dataframe():
    df = pd.DataFrame({"sample": []})
    result = rurils.map_tissue_types(df)

    assert result.empty

def test_non_default_sample_col():
    df = pd.DataFrame({"id": ["X_A_Y", "X_B_Y"], "other_col": [1, 2]})
    result = rurils.map_tissue_types(df, sample_col="id", tissue_col="mapped_tissue")

    assert list(result["mapped_tissue"]) == ["Bud stage 1", "Bud stage 2"]
    assert "id" in result.columns
    assert "other_col" not in result.columns

def test_dataframe_integrity():
    df = pd.DataFrame({"sample": ["X_A_Y"], "keep_this": ["yes"]})
    result = rurils.map_tissue_types(df)

    assert "keep_this" not in result.columns  # Verify that only 'sample' and 'tissue' columns are present
    assert "sample" in result.columns and "tissue" in result.columns  # Verify both necessary columns are present

