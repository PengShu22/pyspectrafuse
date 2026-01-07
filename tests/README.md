# Test Suite for pyspectrafuse

This directory contains the test suite for the pyspectrafuse package.

## Running Tests

### Install test dependencies

```bash
pip install -e ".[test]"
```

### Run all tests

```bash
pytest
```

### Run tests with coverage

```bash
pytest --cov=pyspectrafuse --cov-report=html --cov-report=term-missing
```

### Run specific test file

```bash
pytest tests/test_parquet_utils.py
```

### Run specific test

```bash
pytest tests/test_parquet_utils.py::TestParquetPathHandler::test_init
```

## Test Structure

- `conftest.py`: Shared fixtures and pytest configuration
- `test_parquet_utils.py`: Tests for ParquetPathHandler utility
- `test_msp_utils.py`: Tests for MspUtil utility
- `test_sdrf_utils.py`: Tests for SdrfUtil utility
- `test_cli.py`: Tests for CLI commands
- `test_consensus_strategy.py`: Tests for consensus spectrum strategies

## Continuous Integration

Tests are automatically run on GitHub Actions for:
- Multiple Python versions (3.8, 3.9, 3.10, 3.11)
- Multiple operating systems (Ubuntu, macOS)
- On every push and pull request


