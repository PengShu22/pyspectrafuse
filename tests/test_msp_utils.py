"""Tests for msp_utils module."""
from pathlib import Path
from pyspectrafuse.common.msp_utils import MspUtil


class TestMspUtil:
    """Test cases for MspUtil class."""

    def test_init(self):
        """Test MspUtil initialization."""
        util = MspUtil()
        assert util.namespace is not None

    def test_get_num_peaks(self):
        """Test get_num_peaks method."""
        import numpy as np
        array = np.array([[1, 2], [3, 4], [5, 6]])
        num_peaks = MspUtil.get_num_peaks(array)
        assert num_peaks == 3

    def test_get_msp_dict(self):
        """Test get_msp_dict method."""
        msp_dict = MspUtil.get_msp_dict(
            name="Test",
            mw=1000.5,
            comment="Test comment",
            num_peaks=2,
            mz_arr=[100.0, 200.0],
            intensity_arr=[1000.0, 2000.0]
        )
        assert msp_dict['params']['Name'] == "Test"
        assert msp_dict['params']['MW'] == 1000.5
        assert msp_dict['params']['Num peaks'] == 2
        assert len(msp_dict['mz_array']) == 2
        assert len(msp_dict['intensity_array']) == 2

    def test_get_msp_fmt(self, sample_spectrum_data):
        """Test get_msp_fmt method."""
        import pandas as pd
        row = pd.Series(sample_spectrum_data)
        # Convert arrays to string representation as expected by the function
        row['mz_array'] = str(sample_spectrum_data['mz_array'])
        row['intensity_array'] = str(sample_spectrum_data['intensity_array'])
        
        msp_fmt = MspUtil.get_msp_fmt(row)
        assert "Name:" in msp_fmt
        assert "MW:" in msp_fmt
        assert "Comment:" in msp_fmt
        assert "Num peaks:" in msp_fmt
        assert sample_spectrum_data['peptidoform'] in msp_fmt

    def test_write2msp(self, temp_dir):
        """Test write2msp method."""
        output_path = Path(temp_dir) / "test.msp.txt"
        content = "Name: Test\nMW: 1000.0\nNum peaks: 2\n"
        
        MspUtil.write2msp(str(output_path), content)
        
        assert output_path.exists()
        with open(output_path, 'r') as f:
            written_content = f.read()
        assert content in written_content

    def test_write2msp_append(self, temp_dir):
        """Test write2msp append mode."""
        output_path = Path(temp_dir) / "test_append.msp.txt"
        content1 = "Name: Test1\n"
        content2 = "Name: Test2\n"
        
        MspUtil.write2msp(str(output_path), content1)
        MspUtil.write2msp(str(output_path), content2)
        
        with open(output_path, 'r') as f:
            written_content = f.read()
        assert content1 in written_content
        assert content2 in written_content


