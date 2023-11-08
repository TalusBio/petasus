"""The main CLI entry point."""
import sys

import click
from loguru import logger

from .scripts import (
    scan2rt,
    localize_mods,
    sage2lib,
)


@click.group()
def cli() -> None:
    """Petasus Encapsulates TAlus' Search Utility Scripts

    The subcommands help with various tasks that may be needed before an
    EncyclopeDIA DLIB can be created.
    """
    logger.remove()
    logger.add(
        sys.stderr,
        level="INFO",
        format="{level} {message}",
    )


cli.add_command(scan2rt.scan2rt)
cli.add_command(localize_mods.localize_mods)
cli.add_command(sage2lib.sage2lib)

if __name__ == "__main__":
    cli()
