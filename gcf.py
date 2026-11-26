"""GC mass fraction vs metallicity: shared functions."""
from pathlib import Path

import numpy as np
from astropy.table import Table

DATA = Path.home() / "data" / "gc_fraction"
MW_GC_FITS = Path.home() / "IoA Dropbox/Dr V.A. Belokurov/data/catalogues/gc_catalog_updated.fits"


def bk_energy_boundary(lz):
    """BK in-situ/accreted boundary E_b(Lz) in the E-Lz plane.

    Port of the mask in ~/Dropbox (Personal)/idl/work/globular_clusters/gc_catalog_create.pro.
    Units: lz in 1e3 kpc km/s, E in 1e5 km^2/s^2 (Lz > 0 prograde). In-situ: E < E_b(Lz), |Lz| < 4.4.
    """
    lz = np.asarray(lz, float)
    gse = -1.4 + 0.3 * lz**2          # -0.58 < Lz < 0.58
    disc = -1.325 + 0.075 * lz**2     # Lz > 1
    e_lo, e_hi = -1.4 + 0.3 * 0.58**2, -1.325 + 0.075   # polygon edge 0.58 -> 1
    bridge = e_lo + (e_hi - e_lo) * (lz - 0.58) / (1 - 0.58)
    return np.select([lz <= -0.58, lz < 0.58, lz < 1], [-1.3, gse, bridge], disc)


def bk_in_situ(lz, en):
    """True for in-situ GCs; lz in kpc km/s, en in km^2/s^2 (catalogue units)."""
    lz3, en5 = np.asarray(lz) * 1e-3, np.asarray(en) * 1e-5
    return (en5 < bk_energy_boundary(lz3)) & (np.abs(lz3) < 4.4)


def mw_gcs():
    """MW GC catalogue (Baumgardt masses, Harris [Fe/H]) with BK in-situ flag recomputed."""
    t = Table.read(MW_GC_FITS)
    t.convert_bytestring_to_unicode()
    t["IN_SITU"] = bk_in_situ(t["LZ"], t["EN"])
    return t
