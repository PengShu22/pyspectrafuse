import click
from pyspectrafuse.commands.quantmsio2mgf import quantmsio2mgf
from pyspectrafuse.commands.spectrum2msp import spectrum2msp
from pyspectrafuse import __version__

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
REVISION = "0.1.1"


# Cli returns command line requests
@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def cli():
    """
    pyspectrafuse - Command-line utilities for spectrum clustering and conversion.
    
    This tool provides commands for converting parquet files to MGF format and
    generating MSP format files using various consensus spectrum strategies.
    """


cli.add_command(quantmsio2mgf)
cli.add_command(spectrum2msp)


def main():
    cli()


