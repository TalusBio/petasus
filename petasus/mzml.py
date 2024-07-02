import base64
import zlib

import numpy as np
from lxml import etree


def read(mzml_file: str) -> dict[tuple[np.ndarray, np.ndarray]]:
    """Parse the mzML file

    This parser is much faster than Pyteomics, but doesn't do any error
    checking.

    Parameters
    ----------
    mzml_file : Path-like
        The mzML file to parse.

    Returns
    -------
    Dict of Tuple of (np.ndarray, np.ndarray)
        The mass spectra, indexed by the scan identifier.
    """
    spectra = {}
    for _, elem in etree.iterparse(str(mzml_file), tag="{*}spectrum"):
        spec_id, *arrays = _parse_spectrum(elem)
        spectra[spec_id] = arrays

    return spectra


def _parse_spectrum(elem: etree.Element) -> tuple[str, np.ndarray, np.ndarray]:
    """Parse a mass spectrum.

    Parameters
    ----------
    elem : lxml.etree.Element
        The element with the Spectrum tag.

    Returns
    -------
    spec_id : str
        The spectrum identifier.
    mz_array : np.ndarray
        The m/z values.
    intensity_array : np.ndarray
        The intensity values.
    """
    spec_id = elem.get("id")
    bin_list = next(elem.iter("{*}binaryDataArrayList"))

    # spectrum is m/z array and intensity array.
    spec = [None, None]
    for bin_array in bin_list:
        for child in bin_array:
            if child.tag.endswith("cvParam"):
                if child.get("accession") == "MS:1000514":
                    idx = 0
                elif child.get("accession") == "MS:1000515":
                    idx = 1

            if child.tag.endswith("binary"):
                decoded = zlib.decompress(
                    base64.b64decode(child.text.encode("ascii"))
                )
                array = np.frombuffer(bytearray(decoded), dtype=np.float64)

        spec[idx] = array

    return spec_id, *spec
