"""The main CLI entry point."""
import sys

import click
from loguru import logger

from . import scan2rt
from . import localize_mods


@click.group()
def cli() -> None:
    """Petasus Encapsulates TAlus' Search Utility Scripts

    The subcommands help with various tasks that may be needed before an
    EncyclopeDIA DLIB can be created.
    """
    logger.add(sys.stderr, format="{level} {message}")


cli.add_command(scan2rt.scan2rt)
cli.add_command(localize_mods.localize_mods)
