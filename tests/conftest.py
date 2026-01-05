"""Pytest configuration and shared fixtures."""
import pytest
import tempfile
import shutil
from pathlib import Path
import pandas as pd
import numpy as np


@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing."""
    temp_path = tempfile.mkdtemp()
    yield temp_path
    shutil.rmtree(temp_path)


@pytest.fixture
def sample_parquet_dir(temp_dir):
    """Create a sample parquet directory structure."""
    parquet_dir = Path(temp_dir) / "parquet_files"
    parquet_dir.mkdir(parents=True, exist_ok=True)
    return str(parquet_dir)


@pytest.fixture
def sample_sdrf_file(temp_dir):
    """Create a sample SDRF file."""
    sdrf_path = Path(temp_dir) / "sample.sdrf.tsv"
    sdrf_data = pd.DataFrame({
        'comment[data file]': ['file1.mzML', 'file2.mzML'],
        'Characteristics[organism]': ['Homo sapiens', 'Homo sapiens'],
        'comment[instrument]': ['NT=Orbitrap Fusion Lumos;', 'NT=Orbitrap Fusion Lumos;']
    })
    sdrf_data.to_csv(sdrf_path, sep='\t', index=False)
    return str(sdrf_path)


@pytest.fixture
def sample_spectrum_data():
    """Create sample spectrum data for testing."""
    return {
        'usi': 'mzspec:PXD000001:file1.mzML:scan:1',
        'pepmass': 1000.5,
        'charge': 2,
        'mz_array': [100.0, 200.0, 300.0],
        'intensity_array': [1000.0, 2000.0, 3000.0],
        'peptidoform': 'PEPTIDE',
        'posterior_error_probability': 0.01,
        'Nreps': 5
    }


@pytest.fixture
def sample_cluster_df():
    """Create a sample cluster DataFrame."""
    import numpy as np
    return pd.DataFrame({
        'cluster_accession': ['cluster1', 'cluster1', 'cluster2'],
        'usi': ['mzspec:PXD000001:file1.mzML:scan:1', 'mzspec:PXD000001:file1.mzML:scan:2', 'mzspec:PXD000001:file2.mzML:scan:1'],
        'pepmass': [1000.0, 1000.1, 2000.0],
        'charge': [2, 2, 3],
        'mz_array': [
            np.array([100.0, 200.0]),
            np.array([100.1, 200.1]),
            np.array([300.0, 400.0])
        ],
        'intensity_array': [
            np.array([1000.0, 2000.0]),
            np.array([1100.0, 2100.0]),
            np.array([3000.0, 4000.0])
        ],
        'peptidoform': ['PEPTIDE1', 'PEPTIDE1', 'PEPTIDE2'],
        'posterior_error_probability': [0.01, 0.02, 0.03],
        'global_qvalue': [0.001, 0.002, 0.003]
    })

