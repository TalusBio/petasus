"""Converts a search results to a dlib file"""
from pathlib import Path

import click
import polars as pl
from loguru import logger
from ms2ml import Config, Peptide
from ms2ml.data.adapters.mzml import MZMLAdapter
from ms2ml.data.parsing.encyclopedia import write_encyclopedia


def read_peptides(peptides, qvalue):
    """Read the peptides into a polars DataFrame.

    Parameters
    ----------
    peptides : Pathlike
        The Sage or Mokapot tsv file.
    qvalue : float
        The maximum q-value threshold.

    Returns
    -------
    polars.DataFrame
        The parsed and filtered peptides.
    """
    peptide_df = pl.read_csv(peptides, separator="\t")
    peptide_df.columns = [c.lower() for c in peptide_df.columns]
    if "mokapot q-value" in peptide_df.columns:
        qval_col = "mokapot q-value"
    else:
        qval_col = "peptide_q"

    return (
        peptide_df
        .filter(pl.col(qval_col) <= qvalue)
        .select(["peptide", "charge", "filename", "scannr"])
    )


def yield_psms(peptide_df, mzml_file, config):
    """Iterate over confident annotated spectra.

    Parameters
    ----------
    peptide_df : polars.DataFrame
        The parsed and filtered search results.
    mzml_file : PathLike
        One of the mzML files denoted in the "filename" column of the tsv file.
    config: mz2ml.Config
        The ms2ml configuration.

    Yields
    ------
    AnnotatedPeptideSpectrum
        The annotated mass spectrum to write.
    """
    logger.info(f"Parsing PSMs from {mzml_file}")
    fname = Path(mzml_file).name
    adapter = MZMLAdapter(mzml_file, config)

    # This filter is not super robust, but should be good enough.
    psms = peptide_df.filter(pl.col("filename") == fname)
    for seq, charge, fname, scannr in psms.rows():
        peptide = Peptide.from_proforma_seq(f"{seq}/{charge}", config=config)
        yield adapter[scannr].annotate(peptide)


@click.command()
@click.option(
    "-q",
    "--qvalue",
    help="The q-value threshold for library building.",
    default=0.01,
)

@click.argument("peptides")
@click.argument("mzml_files", nargs=-1)
def sage2lib(peptides, mzml_files, qvalue):
    """Convert search results to a spectral library.

    This script converts search results from Sage with or without Mokapot
    postprocessing into an EncyclopeDIA DLIB or Skyline BLIB spectral library.

    PEPTIDES is either the peptide-level output from Sage either with or
    without Mokapot. Multiple runs should always be combined together
    prior to FDR estimation to produce control FDR in the final library.

    The output is a DLIB or BLIB in the current working directory with
    the same stem as the PEPTIDES file.

    """
    peptide_df = read_peptides(peptides, qvalue)
    config = Config(
        mod_mode="delta_mass",
        mod_fixed_mods=(),
        encoding_mod_alias={},
    )

    for mzml_file in mzml_files:
        write_encyclopedia(
            f"{Path(peptides).stem}.dlib",
            yield_psms(peptide_df, mzml_file, config),
            source_file=mzml_file,
        )

    logger.info("DONE!")
