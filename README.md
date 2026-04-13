# pyspectrafuse

[![Tests](https://github.com/bigbio/pyspectrafuse/actions/workflows/tests.yml/badge.svg)](https://github.com/bigbio/pyspectrafuse/actions/workflows/tests.yml)
[![Containers](https://github.com/bigbio/pyspectrafuse/actions/workflows/pyspectrafuse-containers.yml/badge.svg)](https://github.com/bigbio/pyspectrafuse/actions/workflows/pyspectrafuse-containers.yml)
[![codecov](https://codecov.io/gh/bigbio/pyspectrafuse/branch/master/graph/badge.svg)](https://codecov.io/gh/bigbio/pyspectrafuse)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

Python library and CLI for the [spectrafuse](https://github.com/bigbio/spectrafuse) pipeline -- spectral clustering and consensus library generation for mass spectrometry proteomics data.

## Features

- **Parquet to Dat Conversion**: Convert QPX parquet files directly to MaRaCluster's binary `.dat` format (15x smaller than MGF, ~100 bytes/spectrum)
- **Cluster DB Builder**: Generate `cluster_metadata.parquet` and `psm_cluster_membership.parquet` from MaRaCluster output
- **Consensus Spectrum Generation**: Multiple strategies for generating consensus spectra from clustered data:
  - Best spectrum selection (lowest PEP/q-value)
  - Binning-based consensus
  - Most similar spectrum
  - Average spectrum aggregation
- **Incremental Clustering**: Add new datasets to existing cluster DBs without re-clustering
- **QPX Metadata Handling**: Reads species/instrument from `.run.parquet` + `.sample.parquet` (QPX format)
- **MSP Format Export**: Export consensus spectra to MSP (MassBank) format

## Installation

### Using pip

```bash
pip install pyspectrafuse
```

### Using conda

```bash
conda env create -f environment.yml
conda activate pyspectrafuse
```

### Development Installation

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

#### Convert Parquet to Dat (recommended)

Convert QPX parquet files directly to MaRaCluster's binary `.dat` format:

```bash
pyspectrafuse convert-dat \
    -p /path/to/project/ \
    -o dat_output/ \
    -c 2  # optional: filter by charge state
```

The `.dat` format is 15x smaller than MGF and MaRaCluster reads it directly with the `-D` flag, skipping its own file conversion step.

#### Build Cluster DB

Generate cluster metadata and PSM membership parquets from MaRaCluster output:

```bash
pyspectrafuse build-cluster-db \
    --cluster_tsv clusters_p30.tsv \
    --scan_titles dat_output/project.psm_charge2.scan_titles.txt \
    --parquet_dir /path/to/project/ \
    --dataset_name PXD014877 \
    --species "Homo sapiens" \
    --instrument "Q Exactive HF" \
    --charge charge2 \
    --output_dir cluster_db/
```

#### Generate MSP Consensus Spectra

```bash
pyspectrafuse msp \
    --parquet_dir /path/to/project \
    --method_type best \
    --cluster_tsv_file /path/to/cluster.tsv \
    --species "Homo sapiens" \
    --instrument "Q Exactive HF" \
    --charge "charge2"
```

**Consensus methods:** `best` (default), `bin`, `most`, `average`

#### Incremental Clustering

Add new datasets to an existing cluster DB:

```bash
# Step 1: Extract representative spectra from existing clusters
pyspectrafuse incremental extract-reps \
    --cluster_metadata cluster_db/cluster_metadata.parquet \
    --output_mgf reps.mgf

# Step 2: Merge new data with existing clusters
pyspectrafuse incremental merge-clusters \
    --cluster_tsv new_clusters.tsv \
    --rep_mgf reps.mgf \
    --existing_metadata cluster_db/cluster_metadata.parquet \
    --existing_membership cluster_db/psm_cluster_membership.parquet \
    --new_parquet_dir /path/to/new_project/ \
    --species "Homo sapiens" \
    --instrument "Q Exactive HF" \
    --charge charge2 \
    --method_type bin
```

#### Convert MSNet to QPX

One-time data preparation for legacy MSNet parquet files:

```bash
pyspectrafuse convert-to-qpx \
    --input data.parquet \
    --sdrf data.sdrf.tsv \
    --output data.psm.parquet
```

### Python API

```python
from pyspectrafuse.maracluster_dat import parquet_to_dat
from pyspectrafuse.common.qpx_metadata import get_metadata_dict
from pyspectrafuse.common.parquet_utils import ParquetPathHandler

# Convert parquet to dat (recommended path)
written, skipped, dat_path = parquet_to_dat(
    'data/PXD014877/PXD014877.psm.parquet',
    output_dir='dat_files/',
    file_idx=0,
    charge_filter=2,
)

# Get metadata from QPX parquet files
metadata = get_metadata_dict('/path/to/project/')

# List PSM parquet files
parquet_files = list(ParquetPathHandler.iter_parquet_dir('/path/to/project/'))
```

## Project Structure

```
pyspectrafuse/
├── maracluster_dat.py              # Parquet to .dat binary conversion
├── cluster_parquet_combine/        # Cluster DB builder
├── commands/                       # CLI command implementations
│   ├── build_cluster_db.py         # build-cluster-db command
│   ├── parquet2dat.py              # convert-dat command
│   ├── cluster2parquet.py          # cluster-parquet command
│   ├── quantmsio2mgf.py            # convert-mgf command (legacy)
│   ├── spectrum2msp.py             # msp command
│   ├── incremental.py              # incremental sub-commands
│   └── convert_msnet_to_qpx.py     # convert-to-qpx command
├── common/                         # Shared utilities
│   ├── qpx_metadata.py             # QPX metadata reader
│   ├── parquet_utils.py            # Parquet file handling
│   ├── msp_utils.py                # MSP format utilities
│   └── sdrf_utils.py               # SDRF file parsing (legacy)
├── consensus_strategy/             # Consensus spectrum strategies
│   ├── best_spetrum_strategy.py
│   ├── binning_strategy.py
│   ├── most_similar_strategy.py
│   └── average_spectrum_strategy.py
├── incremental/                    # Incremental clustering modules
│   ├── representative_mgf.py
│   ├── resolve_clusters.py
│   └── merge_results.py
└── mgf_convert/                    # MGF conversion (legacy)
```

## Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=pyspectrafuse --cov-report=html

# Run specific test file
pytest tests/test_parquet_utils.py
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

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
  version = {0.0.4},
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
