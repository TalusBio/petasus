"""Get the retention times, given scan numbers."""
import gzip
import logging
import re
from io import BytesIO

import click
import polars as pl
from loguru import logger
from pyteomics.mzml import MzML

from .. import utils


def parse_mzml(mzml_data):
    """Parse the mzML file

    Parameters
    ----------
    mzml_data : file-like
        The mzML file data.

    Returns
    -------
    Dict
        A dictionary mapping scans to retention times.
    """
    rt_map = []
    with MzML(mzml_data) as run:
        total = len(run)
        logging.info("Reading %i scans...", total)
        for idx, spectrum in enumerate(run):
            scan = (
                re.match(r".*(scan|index)=(.+)$", spectrum["id"])
                .groups()[1]
            )
            rt = spectrum["scanList"]["scan"][0]["scan start time"]
            rt_map.append((int(scan), float(rt)))
            if not (idx + 1) % 10000:
                logger.info("{} / {}", idx + 1, total)

    return pl.DataFrame(rt_map, schema=["scan", "rt"])


@click.command()
@click.argument("mzml_files", nargs=-1)
def scan2rt(mzml_files):
    """Map scan numbers to retention times.

    MZML_FILES should be one or more mzML files, either gzip compressed
    or not.

    This tool outputs a Parquet file mapping scan numbers to retention times.
    """
    for mzml_file in utils.listify(mzml_files):
        # Decompress the file if it is gzipped:
        try:
            with gzip.open(mzml_file, "rb") as gzipped:
                logger.info("Decompressing {}...", mzml_file)
                mzml_data = BytesIO(gzipped.read())
        except gzip.BadGzipFile:
            mzml_data = mzml_file

        map_df = parse_mzml(mzml_data)
        out_file = mzml_file.replace(".gz", "")
        out_file = mzml_file.replace(".mzML", ".scan2rt.parquet")
        logger.info("Writing %s...", out_file)
        map_df.write_parquet(out_file)

    logging.info("DONE!")
