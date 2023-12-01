"""Test the localize mods functionality"""
import math

import numpy as np
from petasus import masses
from petasus.scripts import localize_mods


def test_shift_fragments():
    """Test that fragment calculations are working correctly."""
    frags = np.arange(6).reshape(2, -1).T
    mod = 20.0

    b_shift = localize_mods.shift_fragments(frags, "b", mod)
    b_out = np.array(
        [
            [0.0, 1, 2, 3, 4, 5],
            [20.0, 1, 2, 13, 4, 5],
            [20, 21, 2, 13, 14, 5],
            [20, 21, 22, 13, 14, 15],
        ]
    )

    np.testing.assert_allclose(b_shift, b_out)

    y_shift = localize_mods.shift_fragments(frags, "y", mod)
    y_out = np.array(
        [
            [20.0, 21, 22, 13, 14, 15],
            [0, 21, 22, 3, 14, 15],
            [0, 1, 22, 3, 4, 15],
            [0, 1, 2, 3, 4, 5],
        ]
    )

    np.testing.assert_allclose(y_shift, y_out)


def test_log_factorial():
    """Test that the log factorial function works."""
    n_array = np.arange(5)
    np.testing.assert_allclose(
        localize_mods.log_factorial(n_array),
        np.log([math.factorial(n) for n in n_array]),
    )


def test_shifted_hyperscores():
    """Test that the hyperscores are working.

    Note that the best way to test this is to compute the
    hyperscore manually between a spectrum and peptide,
    but I'm feeling lazy.
    """
    peptide = "LESLIEK"
    mod_peptide = "LES[+79]LIEK"
    mz_array = np.concatenate(
        masses.calculate_fragments(mod_peptide, 2)
    ).flatten()

    mz_array.sort()
    int_array = np.ones(mz_array.shape)
    scores = localize_mods.shifted_hyperscores(
        *masses.calculate_fragments(peptide, 2),
        mz_array=mz_array,
        int_array=int_array,
        mass=79.0,
        tol=(-10, 10),
    )

    assert len(scores) == len(peptide)
    assert np.argmax(scores) == 2
