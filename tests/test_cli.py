"""Tests for CLI commands."""
from click.testing import CliRunner
from pyspectrafuse.pyspectrafuse import cli


class TestCLI:
    """Test cases for CLI commands."""

    def test_cli_help(self):
        """Test CLI help command."""
        runner = CliRunner()
        result = runner.invoke(cli, ['--help'])
        assert result.exit_code == 0
        assert 'convert-dat' in result.output.lower() or 'msp' in result.output.lower()

    def test_cli_version(self):
        """Test CLI version command."""
        runner = CliRunner()
        result = runner.invoke(cli, ['--version'])
        assert result.exit_code == 0
        assert 'version' in result.output.lower()

    def test_convert_dat_command_help(self):
        """Test convert-dat command help."""
        runner = CliRunner()
        result = runner.invoke(cli, ['convert-dat', '--help'])
        assert result.exit_code == 0
        assert 'parquet' in result.output.lower() or 'dat' in result.output.lower()

    def test_msp_command_help(self):
        """Test msp command help."""
        runner = CliRunner()
        result = runner.invoke(cli, ['msp', '--help'])
        assert result.exit_code == 0
        assert 'msp' in result.output.lower() or 'parquet' in result.output.lower()

    def test_convert_dat_command_missing_args(self):
        """Test convert-dat command with missing required arguments."""
        runner = CliRunner()
        result = runner.invoke(cli, ['convert-dat'])
        assert result.exit_code != 0 or '--parquet_dir' in result.output

    def test_msp_command_missing_args(self):
        """Test msp command with missing required arguments."""
        runner = CliRunner()
        result = runner.invoke(cli, ['msp'])
        assert result.exit_code != 0 or '--parquet_dir' in result.output

    def test_incremental_help(self):
        """Test incremental command help."""
        runner = CliRunner()
        result = runner.invoke(cli, ['incremental', '--help'])
        assert result.exit_code == 0
        assert 'extract-reps-dat' in result.output
        assert 'merge-clusters' in result.output

    def test_incremental_extract_reps_dat_help(self):
        """Test incremental extract-reps-dat command help."""
        runner = CliRunner()
        result = runner.invoke(cli, ['incremental', 'extract-reps-dat', '--help'])
        assert result.exit_code == 0
        assert 'cluster_metadata' in result.output

    def test_incremental_merge_clusters_help(self):
        """Test incremental merge-clusters command help."""
        runner = CliRunner()
        result = runner.invoke(cli, ['incremental', 'merge-clusters', '--help'])
        assert result.exit_code == 0
        assert 'scan_titles_dir' in result.output
