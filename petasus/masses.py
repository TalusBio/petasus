"""Utilities to calculate peptide and fragment masses.

Note that this module is particularly tailored to work with Sage's peptide
strings. Your mileage with other search engines may vary!
"""

import re

import numba as nb
import numpy as np

# For parsing Sage peptide strings.
SEQ_REGEX = r"([A-Z]|[+-][\d\.]*)"

# Constants
H = 1.007825035
O = 15.99491463  # noqa: E741
OH = O + H
H2O = 2 * H + O
PROTON = 1.00727646688
C13 = 1.003355
RESIDUES = {
    "G": 57.021463735,
    "A": 71.037113805,
    "S": 87.032028435,
    "P": 97.052763875,
    "V": 99.068413945,
    "T": 101.047678505,
    "C": 103.009184505,  # Note this is without +57.02146
    "L": 113.084064015,
    "I": 113.084064015,
    "N": 114.042927470,
    "D": 115.026943065,
    "Q": 128.058577540,
    "K": 128.094963050,
    "E": 129.042593135,
    "M": 131.040484645,
    "H": 137.058911875,
    "F": 147.068413945,
    "U": 150.953633405,
    "R": 156.101111050,
    "Y": 163.063328575,
    "W": 186.079312980,
    "O": 237.147726925,
}


def calculate_fragments(
    seq: str,
    charge: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Parse a string into a list of masses

    Parameters
    ----------
    seq : str
        The peptide sequence, with modifications. Modification must be
        denoted using a '+' or '-' followed by the modification mass.
        For example, 'LES[+79]IEK' will work.
    charge : int
        The charge of the precursor.

    Returns
    -------
    b_ions
    y_ions : np.ndarray of shape (len(seq) - 1, min(2, charge))
        The b and y ions.
    """
    seq_arr = []
    for aa in re.findall(SEQ_REGEX, seq):
        try:
            seq_arr.append(RESIDUES[aa])
        except KeyError:
            seq_arr[-1] += float(aa)

    return _calc_fragment_masses(np.array(seq_arr), charge)


@nb.njit
def _calc_fragment_masses(
    seq: np.ndarray,
    charge: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate the b and y ions for a peptide sequence.

    Parameters
    ----------
    seq : np.ndarray
        The sequence of residue masses.
    charge : int, optional
        The precursor charge state to consider. If 1, only +1 fragment ions
        are returned. Otherwise, +2 fragment ions are returned.


    Returns
    -------
    b_ions
    y_ions : np.ndarray of shape (len(seq) - 1, min(2, charge))
        The b and y ions.
    """
    max_charge = min(charge, 2)
    n_ions = len(seq) - 1

    b_ions = np.empty((n_ions, max_charge))
    y_ions = np.empty((n_ions, max_charge))

    b_mass = 0
    y_mass = H2O
    for idx in range(n_ions):
        b_mass += seq[idx]
        y_mass += seq[-(idx + 1)]
        for cur_charge in range(1, max_charge + 1):
            z_idx = cur_charge - 1
            b_ions[idx, z_idx] = mass2mz(b_mass, cur_charge)
            y_ions[idx, z_idx] = mass2mz(y_mass, cur_charge)

    return b_ions, y_ions


@nb.njit
def mass2mz(mass: float, charge: int) -> float:
    """Calculate the m/z
    Parameters
    ----------
    mass : float
        The neutral mass.
    charge : int
        The charge.
    Returns
    -------
    float
       The m/z
    """
    return (mass / charge) + PROTON
