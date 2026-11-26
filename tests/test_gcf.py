import numpy as np

import gcf


def test_bk_boundary_matches_idl_flag():
    t = gcf.mw_gcs()
    bad = t["IN_SITU"] != (t["FLAG"] == 1)
    assert bad.sum() == 0, list(t["NAME"][bad])


def test_bk_boundary_continuous():
    lz = np.array([-0.58, 0.58, 1.0])
    left = gcf.bk_energy_boundary(lz - 1e-9)
    right = gcf.bk_energy_boundary(lz + 1e-9)
    assert np.allclose(left, right, atol=1e-3)


def test_mw_names_matchable():
    assert (np.asarray(gcf.mw_gcs()["NAME"]) == "NGC 5139").sum() == 1
