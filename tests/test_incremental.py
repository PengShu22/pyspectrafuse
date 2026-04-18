"""Tests for incremental clustering modules (dat-based)."""
import struct
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from pyspectrafuse.incremental.representative_dat import extract_representatives_dat
from pyspectrafuse.incremental.resolve_clusters_dat import (
    parse_all_scan_titles, resolve_incremental_clusters, parse_cluster_tsv)
from pyspectrafuse.maracluster_dat import SPECTRUM_STRUCT, SCANINFO_STRUCT


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


class TestExtractRepresentativesDat:
    def test_basic_extraction(self, sample_cluster_metadata, tmp_path):
        output_dir = str(tmp_path / 'reps_output')
        written, skipped, dat_path, scaninfo_path, titles_path = \
            extract_representatives_dat(sample_cluster_metadata, output_dir)

        assert written == 3
        assert Path(dat_path).exists()
        assert Path(scaninfo_path).exists()
        assert Path(titles_path).exists()

        # Verify .dat file size (100 bytes per spectrum)
        assert Path(dat_path).stat().st_size == 3 * SPECTRUM_STRUCT.size

        # Verify scan_titles contain rep: prefix
        titles = Path(titles_path).read_text()
        assert 'rep:CLUS_001' in titles
        assert 'rep:CLUS_002' in titles
        assert 'rep:CLUS_003' in titles

    def test_charge_filter(self, sample_cluster_metadata, tmp_path):
        output_dir = str(tmp_path / 'reps_filtered')
        written, skipped, dat_path, _, titles_path = \
            extract_representatives_dat(
                sample_cluster_metadata, output_dir, charge_filter=2)

        # Should only extract charge 2 clusters (CLUS_001, CLUS_002)
        assert written == 2
        titles = Path(titles_path).read_text()
        assert 'rep:CLUS_001' in titles
        assert 'rep:CLUS_002' in titles
        assert 'rep:CLUS_003' not in titles

    def test_max_clusters(self, sample_cluster_metadata, tmp_path):
        output_dir = str(tmp_path / 'reps_limited')
        written, _, _, _, _ = extract_representatives_dat(
            sample_cluster_metadata, output_dir, max_clusters=2)

        assert written == 2

    def test_dat_struct_format(self, sample_cluster_metadata, tmp_path):
        output_dir = str(tmp_path / 'reps_struct')
        written, _, dat_path, _, _ = \
            extract_representatives_dat(sample_cluster_metadata, output_dir)

        # Read and unpack first spectrum
        with open(dat_path, 'rb') as f:
            data = f.read(SPECTRUM_STRUCT.size)
        fields = SPECTRUM_STRUCT.unpack(data)

        assert fields[0] == 0  # file_idx
        assert fields[1] == 0  # scannr (first spectrum)
        assert fields[2] == 2  # charge (CLUS_001 has charge 2)
        assert abs(fields[3] - 500.25) < 0.01  # precursor_mz


class TestParseScanTitles:
    def test_basic_parsing(self, tmp_path):
        titles_content = (
            "0\t0\trep:CLUS_001\n"
            "0\t1\trep:CLUS_002\n"
            "0\t2\tid=mzspec::file.raw:scan:123:PEPTIDE/2\n"
        )
        titles_path = str(tmp_path / 'test.scan_titles.txt')
        Path(titles_path).write_text(titles_content)

        df = parse_all_scan_titles([titles_path])

        assert len(df) == 3
        # Check rep identification
        assert df.iloc[0]['is_rep'] == True
        assert df.iloc[0]['old_cluster_id'] == 'CLUS_001'
        assert df.iloc[2]['is_rep'] == False
        assert df.iloc[2]['old_cluster_id'] is None


class TestResolveIncrementalClustersDat:
    def test_basic_resolution(self, tmp_path):
        # Create representative scan_titles
        rep_titles = (
            "0\t0\trep:OLD_001\n"
            "0\t1\trep:OLD_002\n"
        )
        rep_path = str(tmp_path / 'reps.scan_titles.txt')
        Path(rep_path).write_text(rep_titles)

        # Create new data scan_titles
        new_titles = (
            "0\t0\tid=mzspec::file.raw:scan:1:PEPTIDE/2\n"
            "0\t1\tid=mzspec::file.raw:scan:2:ANOTHER/2\n"
        )
        new_path = str(tmp_path / 'new_data.scan_titles.txt')
        Path(new_path).write_text(new_titles)

        # Create a MaRaCluster TSV where:
        # - rep spectrum 0 → new cluster NEW_A
        # - rep spectrum 1 → new cluster NEW_B
        # - new spectrum 0 → NEW_A (merges with OLD_001)
        # - new spectrum 1 → NEW_C (entirely new cluster)
        tsv_content = (
            "reps.mgf\t0\tNEW_A\n"
            "reps.mgf\t1\tNEW_B\n"
            "new_data.mgf\t0\tNEW_A\n"
            "new_data.mgf\t1\tNEW_C\n"
        )
        tsv_path = str(tmp_path / 'clusters.tsv')
        Path(tsv_path).write_text(tsv_content)

        id_map, new_spectra_df = resolve_incremental_clusters(
            tsv_path, [rep_path, new_path])

        # NEW_A should map to OLD_001
        assert id_map['NEW_A'] == 'OLD_001'
        # NEW_B should map to OLD_002
        assert id_map['NEW_B'] == 'OLD_002'
        # NEW_C has no representative -> UUID
        assert id_map['NEW_C'] != 'NEW_C'
        assert len(id_map['NEW_C']) == 36  # UUID format

        # new_spectra_df should only contain non-representative spectra
        assert len(new_spectra_df) == 2
        assert new_spectra_df.iloc[0]['resolved_cluster_id'] == 'OLD_001'

    def test_cluster_merge(self, tmp_path):
        """Test when two old clusters merge into one new cluster."""
        rep_titles = (
            "0\t0\trep:OLD_001\n"
            "0\t1\trep:OLD_002\n"
        )
        rep_path = str(tmp_path / 'reps.scan_titles.txt')
        Path(rep_path).write_text(rep_titles)

        # Both reps land in same new cluster
        tsv_content = (
            "reps.mgf\t0\tMERGED\n"
            "reps.mgf\t1\tMERGED\n"
        )
        tsv_path = str(tmp_path / 'clusters.tsv')
        Path(tsv_path).write_text(tsv_content)

        id_map, _ = resolve_incremental_clusters(tsv_path, [rep_path])

        # Merged cluster should use the first old ID
        assert id_map['MERGED'] == 'OLD_001'
