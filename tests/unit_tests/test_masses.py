"""Test that mass calculations work correctly."""
import numpy as np
from petasus import masses


# Calculated using Pyteomics:
# These are [b_ions, y_ions]
LESLIEK_PLUS_ONE = [
    [
        114.09134044390001,
        243.13393353187,
        330.16596193614004,
        443.25002591327006,
        556.3340898904,
        685.37668297837,
    ],
    [
        147.11280416447,
        276.15539725243997,
        389.23946122957,
        502.3235252067,
        589.3555536109699,
        718.39814669894,
    ],
]

LESLIEK_PLUS_TWO = [
    [
        57.54930845533501,
        122.07060499932,
        165.58661920145502,
        222.12865119002004,
        278.670683178585,
        343.19197972257,
    ],
    [
        74.06004031561999,
        138.581336859605,
        195.12336884817,
        251.66540083673502,
        295.18141503886994,
        359.702711582855,
    ],
]


def test_plus_one_fragments():
    """Test that the plus one fragments are correct."""
    b_ions, y_ions = masses.calculate_fragments("LESLIEK", 1)
    np.testing.assert_allclose(b_ions[:, 0], np.array(LESLIEK_PLUS_ONE[0]))
    np.testing.assert_allclose(y_ions[:, 0], np.array(LESLIEK_PLUS_ONE[1]))


def test_plus_two_fragments():
    """Test that the plus two fragments are correct."""
    b_ions, y_ions = masses.calculate_fragments("LESLIEK", 2)
    np.testing.assert_allclose(b_ions[:, 1], np.array(LESLIEK_PLUS_TWO[0]))
    np.testing.assert_allclose(y_ions[:, 1], np.array(LESLIEK_PLUS_TWO[1]))


def test_with_mod():
    """Test modification handling."""
    b_ions, y_ions = masses.calculate_fragments("LES[+79]LIEK", 1)
    truth = np.array(LESLIEK_PLUS_ONE)

    # Shift the b and y ions:
    truth[0, 2:] = truth[0, 2:] + 79.0
    truth[1, 4:] = truth[1, 4:] + 79.0

    np.testing.assert_allclose(b_ions[:, 0], truth[0, :])
    np.testing.assert_allclose(y_ions[:, 0], truth[1, :])
