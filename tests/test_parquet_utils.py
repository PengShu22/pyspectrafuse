"""Tests for parquet_utils module."""
from pathlib import Path
from pyspectrafuse.common.parquet_utils import ParquetPathHandler


class TestParquetPathHandler:
    """Test cases for ParquetPathHandler class."""

    def test_init(self, temp_dir):
        """Test ParquetPathHandler initialization."""
        handler = ParquetPathHandler(temp_dir)
        assert handler.parquet_path == temp_dir
        assert isinstance(handler.path_obj, Path)

    def test_get_mgf_filename(self, temp_dir):
        """Test get_mgf_filename method."""
        handler = ParquetPathHandler(temp_dir)
        filename = handler.get_mgf_filename(mgf_file_index=1)
        assert filename.endswith('.mgf')
        assert '1' in filename

    def test_get_item_info(self, temp_dir):
        """Test get_item_info method."""
        # Create a path with a project identifier
        project_path = Path(temp_dir) / "PXD123456"
        project_path.mkdir()
        handler = ParquetPathHandler(str(project_path))
        item_info = handler.get_item_info()
        assert item_info == "PXD123456"

    def test_get_item_info_with_dash(self, temp_dir):
        """Test get_item_info with path containing dashes."""
        project_path = Path(temp_dir) / "PXD123456-something"
        project_path.mkdir()
        handler = ParquetPathHandler(str(project_path))
        item_info = handler.get_item_info()
        assert item_info == "PXD123456"

    def test_iter_parquet_dir_no_files(self, temp_dir):
        """Test iter_parquet_dir with no parquet files."""
        result = ParquetPathHandler.iter_parquet_dir(temp_dir)
        assert isinstance(result, list)
        assert len(result) == 0

    def test_iter_parquet_dir_with_files(self, sample_parquet_dir):
        """Test iter_parquet_dir with parquet files."""
        # Create a parquet file in the correct structure
        parquet_file = Path(sample_parquet_dir) / "test.parquet"
        parquet_file.touch()
        
        # Get parent directory
        parent_dir = Path(sample_parquet_dir).parent
        result = ParquetPathHandler.iter_parquet_dir(str(parent_dir))
        assert len(result) > 0
        assert any('test.parquet' in str(p) for p in result)



