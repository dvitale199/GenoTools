from pathlib import Path

import pandas as pd


def _read_bim_snps(prefix: str) -> set:
    bim = pd.read_csv(f"{prefix}.bim", sep=r"\s+", header=None)
    return set(bim[1])


def test_get_common_snps_matches_legacy(synth_ref_bfile, tmp_path):
    """New port yields the same extracted SNP set + common_snps file as legacy."""
    from genotools.ancestry.preprocessing import get_common_snps as new_fn
    from genotools.utils import get_common_snps as legacy_fn

    # Use the same bfile as both inputs (guarantees a nonempty common set)
    g1, g2 = str(synth_ref_bfile), str(synth_ref_bfile)

    legacy_out = tmp_path / "legacy_common"
    new_out = tmp_path / "new_common"

    legacy_res = legacy_fn(g1, g2, str(legacy_out))
    new_res = new_fn(g1, g2, str(new_out))

    assert set(new_res.keys()) == set(legacy_res.keys())
    assert _read_bim_snps(str(new_out)) == _read_bim_snps(str(legacy_out))
    assert (Path(f"{new_out}.common_snps").read_text()
            == Path(f"{legacy_out}.common_snps").read_text())
