"""Test that localize_mods.py works."""
from pathlib import Path

from click.testing import CliRunner
from petasus.scripts.localize_mods import localize_mods


def test_script(data_path, tmp_path):
    """Test the script."""
    args = [
        str(data_path / "results.sage.parquet"),
        str(data_path / "LQSRPAAPPAPGPGQLTLR.mzML"),
        str(data_path / "config.json"),
    ]

    runner = CliRunner()
    result = runner.invoke(localize_mods, args, catch_exceptions=False)
    assert result.exit_code == 0

    assert Path("results.sage.localized.parquet").exists()
