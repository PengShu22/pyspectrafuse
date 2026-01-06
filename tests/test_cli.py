"""Tests for CLI commands."""
from click.testing import CliRunner
from pyspectrafuse.pyspectrafuse_cli import cli


class TestCLI:
    """Test cases for CLI commands."""

    def test_cli_help(self):
        """Test CLI help command."""
        runner = CliRunner()
        result = runner.invoke(cli, ['--help'])
        assert result.exit_code == 0
        # Check for actual content in the help output
        assert 'convert-mgf' in result.output.lower() or 'msp' in result.output.lower()

    def test_cli_version(self):
        """Test CLI version command."""
        runner = CliRunner()
        result = runner.invoke(cli, ['--version'])
        assert result.exit_code == 0
        assert 'version' in result.output.lower()

    def test_convert_mgf_command_help(self):
        """Test convert-mgf command help."""
        runner = CliRunner()
        result = runner.invoke(cli, ['convert-mgf', '--help'])
        assert result.exit_code == 0
        assert 'parquet' in result.output.lower() or 'mgf' in result.output.lower()

    def test_msp_command_help(self):
        """Test msp command help."""
        runner = CliRunner()
        result = runner.invoke(cli, ['msp', '--help'])
        assert result.exit_code == 0
        assert 'msp' in result.output.lower() or 'parquet' in result.output.lower()

    def test_convert_mgf_command_missing_args(self):
        """Test convert-mgf command with missing required arguments."""
        runner = CliRunner()
        result = runner.invoke(cli, ['convert-mgf'])
        # Should fail or show help when required args are missing
        assert result.exit_code != 0 or '--parquet_dir' in result.output

    def test_msp_command_missing_args(self):
        """Test msp command with missing required arguments."""
        runner = CliRunner()
        result = runner.invoke(cli, ['msp'])
        # Should fail or show help when required args are missing
        assert result.exit_code != 0 or '--parquet_dir' in result.output

