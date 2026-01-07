"""Tests for sdrf_utils module."""
import pytest
from pathlib import Path
from pyspectrafuse.common.sdrf_utils import SdrfUtil


class TestSdrfUtil:
    """Test cases for SdrfUtil class."""

    def test_init(self, temp_dir):
        """Test SdrfUtil initialization."""
        util = SdrfUtil(temp_dir)
        assert util.folder == temp_dir

    def test_get_sdrf_file_path_found(self, temp_dir, sample_sdrf_file):
        """Test get_sdrf_file_path when file exists."""
        sdrf_path = SdrfUtil.get_sdrf_file_path(temp_dir)
        assert Path(sdrf_path).exists()
        assert sdrf_path.endswith('.sdrf.tsv')

    def test_get_sdrf_file_path_not_found(self, temp_dir):
        """Test get_sdrf_file_path when file doesn't exist."""
        empty_dir = Path(temp_dir) / "empty"
        empty_dir.mkdir()
        with pytest.raises(FileNotFoundError):
            SdrfUtil.get_sdrf_file_path(str(empty_dir))

    def test_get_metadata_dict_from_sdrf(self, sample_sdrf_file):
        """Test get_metadata_dict_from_sdrf method."""
        metadata_dict = SdrfUtil.get_metadata_dict_from_sdrf(sample_sdrf_file)
        assert isinstance(metadata_dict, dict)
        assert len(metadata_dict) > 0
        # Check that values are lists (organism_instrument)
        for value in metadata_dict.values():
            assert isinstance(value, list)
            assert len(value) == 2  # organism and instrument


