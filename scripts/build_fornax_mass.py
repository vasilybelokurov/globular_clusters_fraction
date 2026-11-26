"""Fornax stellar mass per [Fe/H] bin (formed and living), LMC-style: selection-weighted giant counts / eta.

1. Parent: Gaia DR3 within 2 deg, parallax consistent with 0, PM chi2 < 9.21 around the Fornax mean, RGB box in
   ((BP-RP)_0, G_0), 16.5 <= G_0 <= 19.5. Background density from r_ell > 90'.
2. Spectroscopic [Fe/H]: Battaglia+2006 FLAMES CaT members (VizieR J/A+A/459/423), matched to the parent (1").
   Weight per star = net parent count / spec count in its (r_ell x G_0) cell.
3. eta = selected giants per Msun formed (Kroupa 0.08-100) from PARSEC (lmc-old-population grid, lmcold.isochrones.selected_count)
   in the same box; ages weighted by lmcold.mass_conversion windows (A reference, B, C alternatives).
   [M/H] = [Fe/H] + log10(0.694 10^a + 0.301), a = 0.3 at [Fe/H] <= -1.8 falling linearly to 0 at -1.0.
   Living / formed from lmcold.surviving_mass.survival_grid with the same weights.
Constants: centre, e, PA, D (McConnachie 2012, 2012AJ....144....4M); PM (McConnachie & Venn 2020, 2020RNAAS...4..229M).

Usage: PYTHONPATH=. python scripts/build_fornax_mass.py -> ~/data/gc_fraction/fornax_mass.csv, plots/fornax_mass.png
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sqlutilpy as sqlutil
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u

import gcf

LMC = gcf.Path(__file__).resolve().parents[2] / "lmc-ancient-disc/lmc-old-population"
sys.path.insert(0, str(LMC / "src"))
from lmcold import isochrones, mass_conversion, surviving_mass  # noqa: E402
from lmcold.mass_samples import FEH_EDGES, FEH_CENTRES  # noqa: E402

RA0, DEC0, ELL, PA, DM = 39.997, -34.449, 0.30, 41.0, 5 * np.log10(147e3) - 5
PM0 = np.array([0.381, -0.358])
EBV = 0.02
A_G, A_BP, A_RP = 2.74 * EBV, 3.37 * EBV, 2.04 * EBV
G_LIM = (16.5, 19.5)
R_EDGES = np.array([0, 10, 20, 30, 45, 60, 90])   # arcmin, elliptical
G_EDGES = np.array([16.5, 18.0, 18.5, 19.0, 19.5])
CACHE = gcf.DATA / "fornax_gaia.parquet"


def r_ell(ra, dec):
    x = (ra - RA0) * np.cos(np.deg2rad(DEC0))
    y = dec - DEC0
    t = np.deg2rad(PA)
    a, b = x * np.sin(t) + y * np.cos(t), -x * np.cos(t) + y * np.sin(t)
    return 60 * np.hypot(a, b / (1 - ELL))


def gaia():
    if not CACHE.exists():
        q = f"""select source_id, ra, dec, phot_g_mean_mag g, phot_bp_mean_mag bp, phot_rp_mean_mag rp, parallax,
                parallax_error, pmra, pmdec, pmra_error, pmdec_error, pmra_pmdec_corr
                from gaia_dr3.gaia_source where q3c_radial_query(ra, dec, {RA0}, {DEC0}, 2.0) and phot_g_mean_mag < 20"""
        pd.DataFrame(sqlutil.get(q, asDict=True)).to_parquet(CACHE)
    return pd.read_parquet(CACHE)


grid = pd.read_parquet(LMC / "data/raw/parsec/gaiaEDR3.parquet")
grid = grid[grid.label <= 8]


def rgb_edge(mh, logage):
    iso = grid[np.isclose(grid.MH, mh, atol=.02) & np.isclose(grid.logAge, logage, atol=.01) & (grid.label == 3)]
    iso = iso.sort_values("Gmag")
    return iso.Gmag.to_numpy() + DM, (iso.G_BPmag - iso.G_RPmag).to_numpy()


BLUE, RED = rgb_edge(-2.19, 10.1), rgb_edge(0.0, 9.4)


def inside(c0, g0):
    lo = np.interp(g0, *BLUE) - 0.15
    hi = np.interp(g0, *RED) + 0.15
    return (g0 >= G_LIM[0]) & (g0 <= G_LIM[1]) & (c0 > lo) & (c0 < hi)


def parent():
    d = gaia()
    d["g0"], d["c0"] = d.g - A_G, (d.bp - A_BP) - (d.rp - A_RP)
    d["r"] = r_ell(d.ra, d.dec)
    dp = d[["pmra", "pmdec"]].to_numpy() - PM0
    s1, s2, rho = d.pmra_error ** 2 + .03 ** 2, d.pmdec_error ** 2 + .03 ** 2, d.pmra_pmdec_corr * d.pmra_error * d.pmdec_error
    chi2 = (dp[:, 0] ** 2 * s2 - 2 * dp[:, 0] * dp[:, 1] * rho + dp[:, 1] ** 2 * s1) / (s1 * s2 - rho ** 2)
    ok = (chi2 < 9.21) & (np.abs(d.parallax - 0.007) < 3 * d.parallax_error) & inside(d.c0, d.g0)
    return d[ok.to_numpy()]


def cell_area():
    """Sky area [arcmin^2] of each r_ell annulus inside the 2-deg query circle, plus the background region r > 90'."""
    rng = np.random.default_rng(1)
    n = 4_000_000
    rad = 2.0 * np.sqrt(rng.uniform(size=n))
    phi = rng.uniform(0, 2 * np.pi, n)
    dec = DEC0 + rad * np.sin(phi)
    ra = RA0 + rad * np.cos(phi) / np.cos(np.deg2rad(DEC0))
    total = np.pi * 120 ** 2
    r = r_ell(ra, dec)
    return np.histogram(r, R_EDGES)[0] / n * total, np.mean(r > R_EDGES[-1]) * total


def spec(p):
    b = Table.read(gcf.DATA / "raw/fornax/J_A+A_459_423_table4.fits")
    b.convert_bytestring_to_unicode()
    b = b[b["Member"] == "Y"].to_pandas()
    c = SkyCoord(b.RAJ2000, b.DEJ2000, unit=(u.hourangle, u.deg))
    idx, sep, _ = c.match_to_catalog_sky(SkyCoord(p.ra.to_numpy() * u.deg, p.dec.to_numpy() * u.deg))
    b["feh"] = b["[Fe/H]"].astype(float)
    m = (sep < 1 * u.arcsec) & np.isfinite(b.feh.to_numpy())
    return pd.concat([p.iloc[idx[m]].reset_index(drop=True), b.loc[m, ["feh"]].reset_index(drop=True)], axis=1), len(b)


def weights(p, s):
    area, bg_area = cell_area()
    np_cell = np.histogram2d(p.r, p.g0, [R_EDGES, G_EDGES])[0]
    bg = p[p.r > R_EDGES[-1]]
    bg_cell = np.histogram(bg.g0, G_EDGES)[0][None, :] / bg_area * area[:, None]
    net = np.clip(np_cell - bg_cell, 0, None)
    ns = np.histogram2d(s.r, s.g0, [R_EDGES, G_EDGES])[0]
    w = np.divide(net, ns, out=np.zeros_like(net), where=ns > 0)
    # cells without spec stars: rescale the covered cells of the same r_ell bin (MDF assumed independent of G_0)
    covered = np.where(ns > 0, net, 0).sum(1)
    w *= np.divide(net.sum(1), covered, out=np.zeros_like(covered), where=covered > 0)[:, None]
    i = np.clip(np.digitize(s.r, R_EDGES) - 1, 0, len(R_EDGES) - 2)
    j = np.clip(np.digitize(s.g0, G_EDGES) - 1, 0, len(G_EDGES) - 2)
    missing = net[ns == 0].sum() / net.sum()
    return w[i, j], net, ns, bg_cell, missing


def mh_of(feh):
    a = np.clip(0.3 * (-1.0 - feh) / 0.8, 0, 0.3)
    return feh + np.log10(0.694 * 10 ** a + 0.301)


def eta_and_living(prior):
    rows = []
    for (mh, la), iso in grid.groupby(["MH", "logAge"]):
        iso = iso.sort_values("Mini", kind="stable")
        n = isochrones.selected_count(iso.int_IMF.to_numpy(), (iso.G_BPmag - iso.G_RPmag).to_numpy(),
                                      iso.Gmag.to_numpy() + DM, inside) / isochrones.PARSEC_IMF_NORMALIZATION
        rows.append((mh, la, n))
    eg = pd.DataFrame(rows, columns=["MH", "logAge", "eta"])
    sv = surviving_mass.survival_grid(grid)
    eg = eg.merge(sv[["mh", "logage", "living_fraction"]], left_on=["MH", "logAge"], right_on=["mh", "logage"])
    mhs = np.sort(eg.MH.unique())
    out = []
    for feh in FEH_CENTRES:
        mh = mhs[np.argmin(np.abs(mhs - mh_of(feh)))]
        e = eg[eg.MH == mh].sort_values("logAge")
        w = mass_conversion.age_weights(10 ** (e.logAge.to_numpy() - 9), feh, prior)
        out.append((np.sum(w * e.eta) / w.sum(), np.sum(w * e.eta * e.living_fraction) / np.sum(w * e.eta)))
    return np.array(out)


p = parent()
s, n_b06 = spec(p)
s["w"], net, ns, bg_cell, missing = weights(p, s)
counts = np.histogram(np.clip(s.feh, -9, 9), FEH_EDGES, weights=s.w)[0]
raw = np.histogram(np.clip(s.feh, -9, 9), FEH_EDGES)[0]
res = {}
for prior in "ABC":
    el = eta_and_living(prior)
    res[prior] = (counts / el[:, 0], counts / el[:, 0] * el[:, 1])
# bootstrap over spec stars (weights fixed): statistical error per bin
rng = np.random.default_rng(2)
el_a = eta_and_living("A")
boot = np.array([np.histogram(np.clip(s.feh, -9, 9), FEH_EDGES, weights=s.w)[0] if False else
                 np.histogram(np.clip(s.feh.to_numpy()[k], -9, 9), FEH_EDGES, weights=s.w.to_numpy()[k])[0]
                 for k in (rng.integers(0, len(s), len(s)) for _ in range(500))]) / el_a[:, 0]
tab = pd.DataFrame(dict(feh=FEH_CENTRES, feh_lo=FEH_EDGES[:-1], feh_hi=FEH_EDGES[1:], n_spec=raw, n_weighted=counts,
                        eta_A=el_a[:, 0], formed_mass=res["A"][0], living_mass=res["A"][1],
                        formed_p16=np.percentile(boot, 16, 0), formed_p84=np.percentile(boot, 84, 0),
                        formed_min=np.min([res[k][0] for k in "ABC"], 0), formed_max=np.max([res[k][0] for k in "ABC"], 0),
                        living_fraction=el_a[:, 1]))
tab.to_csv(gcf.DATA / "fornax_mass.csv", index=False)

print(f"parent {len(p)} (net Fornax {net.sum():.0f}, background {bg_cell.sum():.0f}); B06 members {n_b06}, matched {len(s)}")
print("spec per cell:\n", ns.astype(int), "\nnet parent per cell:\n", net.round(0).astype(int))
print(f"net parent in cells without spec: {missing:.3f}")
print(tab.to_string(float_format="%.3g"))
for k in "ABC":
    print(f"prior {k}: formed {res[k][0].sum():.3e}, living {res[k][1].sum():.3e}")

fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
ax[0].scatter(p.c0, p.g0, s=1, c="0.6", label="parent")
ax[0].scatter(s.c0, s.g0, s=4, c=s.feh, cmap="viridis", vmin=-2.5, vmax=-0.5, label="B06 matched")
gg = np.linspace(*G_LIM, 50)
ax[0].plot(np.interp(gg, *BLUE) - .15, gg, "r-", np.interp(gg, *RED) + .15, gg, "r-")
ax[0].set(xlabel="(BP-RP)$_0$", ylabel="G$_0$", xlim=(0.5, 2.5), ylim=(19.7, 16.3))
ax[1].scatter(p.ra, p.dec, s=1, c="0.6")
ax[1].scatter(s.ra, s.dec, s=4, c="tab:red")
ax[1].set(xlabel="RA [deg]", ylabel="Dec [deg]")
ax[2].step(FEH_CENTRES, tab.formed_mass, where="mid", label="formed (A)")
ax[2].step(FEH_CENTRES, tab.living_mass, where="mid", label="living (A)")
ax[2].fill_between(FEH_CENTRES, tab.formed_min, tab.formed_max, step="mid", alpha=.3, label="formed, priors A–C")
ax[2].errorbar(FEH_CENTRES, tab.formed_mass, [tab.formed_mass - tab.formed_p16, tab.formed_p84 - tab.formed_mass], fmt="none", c="k")
ax[2].set(yscale="log", xlabel="[Fe/H] (B06 CaT)", ylabel=r"M [M$_\odot$] per bin", ylim=(1e4, 1e8))
ax[2].legend(fontsize=8)
fig.suptitle(f"Fornax: formed {res['A'][0].sum():.2e}, living {res['A'][1].sum():.2e} M$_\\odot$ (G$_0$ ≤ 19.5 giants)")
fig.tight_layout()
fig.savefig("plots/fornax_mass.png", dpi=120)
