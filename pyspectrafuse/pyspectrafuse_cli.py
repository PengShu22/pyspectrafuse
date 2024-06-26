import click
from pyspectrafuse.commands.quantmsio2mgf import quantmsio2mgf
from pyspectrafuse.commands.spectrum2msp import spectrum2msp

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
REVISION = "0.1.1"


# Cli returns command line requests
@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version='0.0.24')
def cli():
    """
    xzx
    """


cli.add_command(quantmsio2mgf)
cli.add_command(spectrum2msp)


def main():
    cli()


if __name__ == "__main__":
    main()
