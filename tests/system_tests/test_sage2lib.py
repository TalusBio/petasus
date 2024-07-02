"""Test building a dlib"""
import sqlite3
import subprocess

import pytest


@pytest.mark.parametrize("ftype", ["parquet", "csv"])
def test_simpole_speclib_generation(data_path, tmp_path, ftype):
    """Test spectral library generation"""
    psms = data_path / f"results.sage.{ftype}"
    mzml = data_path / "LQSRPAAPPAPGPGQLTLR.mzML"
    outfile = tmp_path / "results.sage.dlib"
    cmd = [
        "petasus",
        "sage2lib",
        "--qvalue",
        "1",
        f"{psms}",
        f"{mzml}",
    ]

    subprocess.run(cmd, check=True)

    con = sqlite3.Connection(str(outfile))
    cur = con.cursor()
    res = cur.execute("select PeptideModSeq from entries;")
    res = [x[0] for x in res.fetchall()]
    assert len(res)
