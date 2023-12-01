"""Localize Sage modifications.

This script sanitizes the Sage modification masses and attempts to
localize them on the peptide using a crude, highest-scoring method.
"""
import json
import logging
import time
from pathlib import Path

import click
import numba as nb
import numpy as np
import polars as pl
from loguru import logger
from pyteomics.mzml import MzML

from .. import masses


@nb.njit
def shift_fragments(ions, ion_type, mass):
    """Create all possible mass shifts.

    This function abuses numpy broadcasing, so its a bit confusing.

    Parameters
    ----------
    ions : np.ndarray of shape (n_fragments, charge)
        The fragment m/z values, where each column is an increasing
        charge state.
    ion_type : str, {"b", "y"}
        The ion type to shift.
    mass : float
        The mass shift to use.

    Returns
    -------
    np.ndarray of shape (n_fragments, n_fragments * charge)
        The fragment ions where each row are the ions shifted by
        the modification in that position.
    """
    shifts = np.zeros((ions.shape[1], 1, 1))
    shifts[:, 0, 0] = mass / np.arange(1, ions.shape[1] + 1)

    shift_mat = np.zeros((ions.shape[0] + 1, ions.shape[0]))
    shift_mat[1:, :] = np.tri(ions.shape[0])
    if ion_type == "y":
        shift_mat = np.flip(shift_mat)

    shift_mat = shifts * shift_mat

    new_masses = (
        (np.expand_dims(ions.T, 1) + shift_mat)
        .transpose(1, 0, 2)
        .copy()
        .reshape(ions.shape[0] + 1, -1)
    )
    return new_masses


@nb.vectorize([nb.float64(nb.int64)])
def log_factorial(n: int) -> float:
    """Compute the log factorial.

    The major benefit of this function is that the factorial is computed as a
    sum in log space, preventing integer overflows and such.

    Parameters
    ----------
    n : int
        The input integer.

    Returns
    -------
    float
        The log of the n!.
    """
    return np.sum(np.log(np.arange(1, n + 1)))


@nb.njit
def shifted_hyperscores(
    b_ions: np.ndarray,
    y_ions: np.ndarray,
    mz_array: np.ndarray,
    int_array: np.ndarray,
    mass: float,
    tol: tuple[float, float],
) -> np.ndarray:
    """Score every position for the mass shift in a peptide.

    This function is reasonably fast, ~3ms per peptide on my Mac.

    Parameters
    ----------
    b_ions : np.ndarray of shape (n_residues, n_charge_states)
        The m/z values of the b-ion series for a peptide.
    y_ions : np.ndarray of shape (n_residues, n_charge_states)
        The m/z values of the y-ion series for a peptide.
    mz_array : np.ndarray of shape (n_peaks,)
        The m/z values of the mass spectrum.
    int_array : np.ndarray of shape (n_peaks,)
        The intensity values of the mass spectrum.
    mass : float
        The mass shift from the open modification search.
    tol : tuple of float
        The fragment matching tolerances in ppm.

    Returns
    -------
    np.ndarry of shape (n_residues,)
        The log hyperscore the peptide bearing the observed mass shift at
        each amino acid position.
    """
    min_mz = mz_array + (mz_array * tol[0]) / 1e6
    max_mz = mz_array + (mz_array * tol[1]) / 1e6
    ions = np.concatenate(
        (
            shift_fragments(b_ions, "b", mass),
            shift_fragments(y_ions, "y", mass),
        ),
        axis=1,
    )
    ions = np.expand_dims(ions, 2)
    in_tol = (ions > min_mz) & (ions < max_mz)
    matches = np.zeros(in_tol.shape)
    int_array = (
        np.sqrt(int_array)
        .repeat(np.prod(np.array(ions.shape[:2])))
        .reshape(-1, *ions.shape[:2])
        .transpose(1, 2, 0)
    )

    # If Numba supported 3D logical indexing, I would use the following
    # instead of the for loop below:
    # matches[in_tol] = int_array[in_tol]
    for loc in zip(*np.nonzero(in_tol)):
        matches[loc] = int_array[loc]

    # If Numba supported the axis argument in max(), I would use the
    # following instead of the nested for loop below:
    # peaks = matches.max(axis=2)
    peaks = np.empty(matches.shape[:2])
    for row in range(matches.shape[0]):
        for col in range(matches.shape[1]):
            peaks[row, col] = matches[row, col].max()

    # Calculate the Hyperscore in log space:
    num_b = np.sum(peaks[:, : int(peaks.shape[1] / 2)] > 0, axis=1)
    num_y = np.sum(peaks[:, int(peaks.shape[1] / 2) :] > 0, axis=1)
    log_dotp = peaks.sum(axis=1)
    log_dotp[log_dotp > 0] = np.log(log_dotp[log_dotp > 0])
    score = log_dotp + log_factorial(num_b) + log_factorial(num_y)

    # If haven't yet figured out why, but the scores are in
    # reverse order...
    return np.flip(score)


def read_mzml(mzml_file: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Parse an mzML file, extracting the spectra and scan names.

    Parameters
    ----------
    mzml_file : str
        The mzML file to parse.

    Returns
    -------
    Dict of Tuple of (numpy.ndarry, numpy.ndarray)
        The mapping of scan names to mass spectra.
    """
    with MzML(mzml_file) as mzml_data:
        spectra = {
            s["id"]: (s["m/z array"], s["intensity array"]) for s in mzml_data
        }

    return spectra


@click.command("localize-mods")
@click.argument("pin_file")
@click.argument("mzml_file")
@click.argument("config_file")
def localize_mods(pin_file, mzml_file, config_file):
    """Localize a modification from an open modification search.

    We take a simplistic approach: for every position in the peptide
    calculate hyperscore with the modification. The winning position is the one
    that results in the highest hyperscore. Notably, we assume the full
    modification mass may be attributed to a single modification.

    PIN_FILE : A tab-separated or parque file with scannr, expmass, calcmass,
    charge, and peptide columns.

    MZML_FILE : An mzML file, either gzipped or not.

    CONFIG_FILE : The Sage JSON configuration file.

    This tool outputs a modified version of the input PIN file.
    """
    start = time.time()
    logger.info("Reading {}...", pin_file)
    try:
        is_parquet = True
        pin_df = pl.read_parquet(pin_file)
    except pl.ArrowError:
        is_parquet = False
        pin_df = pl.read_csv(
            pin_file,
            separator="\t",
            truncate_ragged_lines=True,
        )

    pin_df = pin_df.with_columns(
        (pl.col("expmass") - pl.col("calcmass")).alias("delta_mass"),
    )

    logger.info("Reading {}...", mzml_file)
    spectra = read_mzml(mzml_file)

    logging.info("Reading Sage configuration...")
    with open(config_file) as conf_json:
        conf = json.load(conf_json)

    tol = tuple(conf["fragment_tol"]["ppm"])

    logger.info("Localizing modifications on {} PSMs...", len(pin_df))
    positions = []
    mod_scores = []
    delta_mod_scores = []
    for idx, psm in enumerate(pin_df.iter_rows(named=True)):
        mz_array, int_array = spectra[psm["scannr"]]
        b_ions, y_ions = masses.calculate_fragments(
            psm["peptide"], psm["charge"]
        )
        scores = shifted_hyperscores(
            b_ions=b_ions,
            y_ions=y_ions,
            mz_array=mz_array,
            int_array=int_array,
            mass=psm["delta_mass"],
            tol=tol,
        )
        top = np.argpartition(-scores, 2)
        positions.append(top[0])
        mod_scores.append(scores[top[0]])
        delta_mod_scores.append(scores[top[0]] - scores[top[1]])
        if not (idx + 1) % 10000:
            percent = int((idx + 1) / len(pin_df) * 100)
            logger.info(" {:0.0f}%", percent)

    logging.info("Writing new PIN file...")
    pin_df = pin_df.with_columns(
        pl.Series(positions).alias("mod_position"),
        pl.Series(mod_scores).alias("shifted_hyperscore"),
        pl.Series(delta_mod_scores).alias("delta_shifted_hyperscore"),
    )

    out_base = Path(pin_file).stem
    if is_parquet:
        pin_df.write_parquet(out_base + ".localized.parquet")
    else:
        pin_df.write_csv(out_base + ".localized.pin", separator="\t")

    done = time.time()
    logger.info("DONE!")
    logger.info(f"Completed in {(done - start) / 60} min.")
