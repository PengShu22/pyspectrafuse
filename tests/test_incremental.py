"""Tests for incremental clustering modules."""
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from pyspectrafuse.incremental.representative_mgf import extract_representatives
from pyspectrafuse.incremental.resolve_clusters import (
    build_representative_index, resolve_incremental_clusters, parse_cluster_tsv)


@pytest.fixture
def sample_cluster_metadata(tmp_path):
    """Create a minimal cluster_metadata.parquet for testing."""
    df = pd.DataFrame({
        'cluster_id': ['CLUS_001', 'CLUS_002', 'CLUS_003'],
        'precursor_mz': [500.25, 600.30, 700.35],
        'charge': [2, 2, 3],
        'consensus_mz_array': [
            [100.0, 200.0, 300.0],
            [150.0, 250.0],
            [180.0, 280.0, 380.0, 480.0],
        ],
        'consensus_intensity_array': [
            [1000.0, 2000.0, 3000.0],
            [1500.0, 2500.0],
            [800.0, 1600.0, 2400.0, 3200.0],
        ],
        'species': ['Homo sapiens'] * 3,
        'instrument': ['Q Exactive'] * 3,
        'peptidoform': ['PEPTIDE/2', 'ANOTHER/2', 'THIRD/3'],
        'peptide_sequence': ['PEPTIDE', 'ANOTHER', 'THIRD'],
        'consensus_method': ['best'] * 3,
        'member_count': [5, 3, 1],
        'project_count': [1, 1, 1],
        'best_pep': [0.001, 0.005, 0.01],
        'best_qvalue': [0.0001, 0.0005, 0.001],
        'purity': [1.0, 0.95, 1.0],
    })

    schema = pa.schema([
        ('cluster_id', pa.string()),
        ('precursor_mz', pa.float64()),
        ('charge', pa.int8()),
        ('consensus_mz_array', pa.list_(pa.float32())),
        ('consensus_intensity_array', pa.list_(pa.float32())),
        ('species', pa.string()),
        ('instrument', pa.string()),
        ('peptidoform', pa.string()),
        ('peptide_sequence', pa.string()),
        ('consensus_method', pa.string()),
        ('member_count', pa.int32()),
        ('project_count', pa.int16()),
        ('best_pep', pa.float64()),
        ('best_qvalue', pa.float64()),
        ('purity', pa.float32()),
    ])

    path = str(tmp_path / 'cluster_metadata.parquet')
    table = pa.Table.from_pandas(df, schema=schema, preserve_index=False)
    pq.write_table(table, path)
    return path


class TestExtractRepresentatives:
    def test_basic_extraction(self, sample_cluster_metadata, tmp_path):
        output = str(tmp_path / 'reps.mgf')
        count = extract_representatives(sample_cluster_metadata, output)

        assert count == 3
        content = Path(output).read_text()
        assert content.count('BEGIN IONS') == 3
        assert 'TITLE=rep:CLUS_001' in content
        assert 'TITLE=rep:CLUS_002' in content
        assert 'TITLE=rep:CLUS_003' in content

    def test_max_clusters(self, sample_cluster_metadata, tmp_path):
        output = str(tmp_path / 'reps_limited.mgf')
        count = extract_representatives(sample_cluster_metadata, output, max_clusters=2)

        assert count == 2
        content = Path(output).read_text()
        assert content.count('BEGIN IONS') == 2

    def test_mgf_format(self, sample_cluster_metadata, tmp_path):
        output = str(tmp_path / 'reps.mgf')
        extract_representatives(sample_cluster_metadata, output)

        content = Path(output).read_text()
        # Check that each record has required MGF fields
        assert 'PEPMASS=' in content
        assert 'CHARGE=' in content
        assert 'END IONS' in content


class TestBuildRepresentativeIndex:
    def test_basic_index(self, tmp_path):
        mgf_content = (
            "BEGIN IONS\n"
            "TITLE=rep:CLUS_001\n"
            "PEPMASS=500.25\n"
            "CHARGE=2+\n"
            "100.0 1000.0\n"
            "END IONS\n"
            "BEGIN IONS\n"
            "TITLE=rep:CLUS_002\n"
            "PEPMASS=600.30\n"
            "CHARGE=2+\n"
            "150.0 1500.0\n"
            "END IONS\n"
        )
        mgf_path = str(tmp_path / 'reps.mgf')
        Path(mgf_path).write_text(mgf_content)

        index = build_representative_index(mgf_path)
        assert index[0] == 'CLUS_001'
        assert index[1] == 'CLUS_002'
        assert len(index) == 2


class TestResolveIncrementalClusters:
    def test_basic_resolution(self, tmp_path):
        # Create a representative MGF with 2 clusters
        mgf_content = (
            "BEGIN IONS\n"
            "TITLE=rep:OLD_001\n"
            "PEPMASS=500.0\n"
            "CHARGE=2+\n"
            "100.0 1000.0\n"
            "END IONS\n"
            "BEGIN IONS\n"
            "TITLE=rep:OLD_002\n"
            "PEPMASS=600.0\n"
            "CHARGE=2+\n"
            "200.0 2000.0\n"
            "END IONS\n"
        )
        mgf_path = str(tmp_path / 'reps.mgf')
        Path(mgf_path).write_text(mgf_content)

        # Create a MaRaCluster TSV where:
        # - rep spectrum 0 → new cluster NEW_A
        # - rep spectrum 1 → new cluster NEW_B
        # - new spectrum in new_data.mgf index 0 → NEW_A (merges with OLD_001)
        # - new spectrum in new_data.mgf index 1 → NEW_C (entirely new cluster)
        tsv_content = (
            f"reps.mgf\t0\tNEW_A\n"
            f"reps.mgf\t1\tNEW_B\n"
            f"new_data.mgf\t0\tNEW_A\n"
            f"new_data.mgf\t1\tNEW_C\n"
        )
        tsv_path = str(tmp_path / 'clusters.tsv')
        Path(tsv_path).write_text(tsv_content)

        id_map, new_spectra_df = resolve_incremental_clusters(tsv_path, mgf_path)

        # NEW_A should map to OLD_001 (representative was in that cluster)
        assert id_map['NEW_A'] == 'OLD_001'
        # NEW_B should map to OLD_002
        assert id_map['NEW_B'] == 'OLD_002'
        # NEW_C has no representative → gets a UUID
        assert id_map['NEW_C'] != 'NEW_C'
        assert len(id_map['NEW_C']) == 36  # UUID format

        # new_spectra_df should only contain non-representative spectra
        assert len(new_spectra_df) == 2
        assert new_spectra_df.iloc[0]['resolved_cluster_id'] == 'OLD_001'

    def test_cluster_merge(self, tmp_path):
        """Test when two old clusters merge into one new cluster."""
        mgf_content = (
            "BEGIN IONS\n"
            "TITLE=rep:OLD_001\n"
            "PEPMASS=500.0\n"
            "CHARGE=2+\n"
            "100.0 1000.0\n"
            "END IONS\n"
            "BEGIN IONS\n"
            "TITLE=rep:OLD_002\n"
            "PEPMASS=500.1\n"
            "CHARGE=2+\n"
            "100.0 1000.0\n"
            "END IONS\n"
        )
        mgf_path = str(tmp_path / 'reps.mgf')
        Path(mgf_path).write_text(mgf_content)

        # Both reps land in same new cluster
        tsv_content = (
            f"reps.mgf\t0\tMERGED\n"
            f"reps.mgf\t1\tMERGED\n"
        )
        tsv_path = str(tmp_path / 'clusters.tsv')
        Path(tsv_path).write_text(tsv_content)

        id_map, _ = resolve_incremental_clusters(tsv_path, mgf_path)

        # Merged cluster should use the first old ID
        assert id_map['MERGED'] == 'OLD_001'
