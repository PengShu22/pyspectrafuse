# pyspectrafuse

Python library with utility scripts for spectrafuse pipeline - a tool for processing and converting mass spectrometry data formats, generating consensus spectra, and handling cluster-based analysis workflows.

## Features

- **Parquet to MGF Conversion**: Convert parquet files to MGF format with parallel processing support
- **Consensus Spectrum Generation**: Multiple strategies for generating consensus spectra from clustered data:
  - Best spectrum selection
  - Most similar spectrum
  - Binning-based consensus
  - Average spectrum aggregation
- **MSP Format Export**: Export spectra to MSP (MassBank) format
- **SDRF File Handling**: Automatic detection and parsing of SDRF (Sample and Data Relationship Format) files
- **Cluster Analysis**: Integration with MaRaCluster output for spectrum clustering workflows

## Installation

### Using pip

```bash
pip install pyspectrafuse
```

### Using conda

```bash
# Create environment from environment.yml
conda env create -f environment.yml
conda activate pyspectrafuse

# Or install directly
conda install -c conda-forge pyspectrafuse
```

### Development Installation

For development, install in editable mode with test dependencies:

```bash
git clone https://github.com/bigbio/pyspectrafuse.git
cd pyspectrafuse
pip install -e ".[test]"
```

## Requirements

- Python >= 3.8
- Click >= 8.0.0
- numpy >= 1.20.0
- pandas >= 1.3.0
- pyarrow >= 5.0.0
- pyteomics >= 4.5.0
- spectrum_utils >= 0.4.0

## Usage

### Command Line Interface

The package provides a command-line interface with two main commands:

#### Convert Parquet to MGF

Convert parquet files to MGF format:

```bash
pyspectrafuse convert-mgf \
    --parquet_dir /path/to/parquet/files \
    --batch_size 100000 \
    --spectra_capacity 1000000 \
    --task_parallel 4
```

**Options:**
- `--parquet_dir, -p`: Directory containing parquet files
- `--batch_size, -b`: Batch size for each parquet pass (default: 100000)
- `--spectra_capacity, -c`: Number of spectra per MGF file (default: 1000000)
- `--task_parallel, -t`: Number of parallel conversion tasks (default: 1)

#### Generate MSP Format Files

Generate MSP format files from clustered spectra:

```bash
pyspectrafuse msp \
    --parquet_dir /path/to/project \
    --method_type average \
    --cluster_tsv_file /path/to/cluster.tsv \
    --species "Homo sapiens" \
    --instrument "Orbitrap Fusion Lumos" \
    --charge "charge2"
```

**Consensus Methods:**
- `best`: Select best spectrum based on posterior error probability
- `most`: Most similar spectrum using similarity measures
- `bin`: Binning-based consensus spectrum
- `average`: Average spectrum aggregation

**Method-specific Options:**

For `most` method:
- `--sim`: Similarity measure (default: 'dot')
- `--fragment_mz_tolerance`: Fragment m/z tolerance (default: 0.02)

For `bin` method:
- `--min_mz`: Minimum m/z (default: 100)
- `--max_mz`: Maximum m/z (default: 2000)
- `--bin_size`: Bin size in m/z (default: 0.02)
- `--peak_quorum`: Peak quorum threshold (default: 0.25)
- `--edge_case_threshold`: Edge case threshold (default: 0.5)

For `average` method:
- `--diff_thresh`: Minimum distance between MS/MS peak clusters (default: 0.01)
- `--dyn_range`: Dynamic range (default: 1000)
- `--min_fraction`: Minimum fraction of cluster spectra (default: 0.5)
- `--pepmass`: Precursor mass calculation method (choices: 'naive_average', 'neutral_average', 'lower_median', default: 'lower_median')
- `--msms_avg`: MS/MS averaging method (choices: 'naive', 'weighted', default: 'weighted')

### Python API

```python
from pyspectrafuse.common.parquet_utils import ParquetPathHandler
from pyspectrafuse.common.sdrf_utils import SdrfUtil
from pyspectrafuse.common.msp_utils import MspUtil

# Get parquet files
parquet_files = ParquetPathHandler.iter_parquet_dir("/path/to/parquet")

# Get SDRF file path
sdrf_path = SdrfUtil.get_sdrf_file_path("/path/to/project")

# Get metadata from SDRF
metadata = SdrfUtil.get_metadata_dict_from_sdrf(sdrf_path)
```

## Testing

Run the test suite:

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=pyspectrafuse --cov-report=html

# Run specific test file
pytest tests/test_parquet_utils.py
```

Test coverage is currently at ~45% and includes tests for:
- CLI commands
- Utility classes (ParquetPathHandler, MspUtil, SdrfUtil)
- Consensus spectrum strategies
- Core functionality

## Project Structure

```
pyspectrafuse/
├── cluster_parquet_combine/    # Cluster and parquet combination utilities
├── commands/                   # CLI command implementations
│   ├── quantmsio2mgf.py       # Parquet to MGF conversion
│   └── spectrum2msp.py        # MSP format generation
├── common/                     # Common utilities
│   ├── msp_utils.py           # MSP format utilities
│   ├── parquet_utils.py        # Parquet file handling
│   └── sdrf_utils.py          # SDRF file parsing
├── consensus_strategy/         # Consensus spectrum strategies
│   ├── average_spectrum_strategy.py
│   ├── best_spetrum_strategy.py
│   ├── binning_strategy.py
│   └── most_similar_strategy.py
└── mgf_convert/                # MGF conversion utilities
```

## Continuous Integration

The project uses GitHub Actions for continuous integration:
- Tests run on Python 3.8, 3.9, 3.10, and 3.11
- Tests run on Ubuntu and macOS
- Coverage reports are generated automatically

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## Citation

If you use pyspectrafuse in your research, please cite:

```bibtex
@software{pyspectrafuse,
  title = {pyspectrafuse: Python tools for spectrafuse pipeline},
  author = {BigBio Team},
  url = {https://github.com/bigbio/pyspectrafuse},
  version = {0.0.2},
  year = {2024}
}
```

## Links

- **Homepage**: https://github.com/bigbio/pyspectrafuse
- **Repository**: https://github.com/bigbio/pyspectrafuse
- **Issues**: https://github.com/bigbio/pyspectrafuse/issues

## Authors

- **BigBio Team** - ypriverol@gmail.com

## Acknowledgments

This project is part of the BigBio initiative for open-source bioinformatics tools.
