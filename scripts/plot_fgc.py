"""f_GC([Fe/H]) = M_GC / M_* for the MW (in-situ GCs, BK boundary), the LMC and Fornax.

GC masses: Baumgardt N-body; MW main curve uses initial masses (MINI), present shown for reference; LMC/Fornax present. Field masses: living stars (MW: mw_mdf k1elfc,
LMC: lmc-old-population giant_present_mass.csv, Fornax: scripts/build_fornax_mass.py); LMC present (living + remnants) shown dashed.
MW field [M/H] (Andrae+23, APOGEE-calibrated) treated as [Fe/H]. GCs binned on each galaxy's field bins.

Usage: PYTHONPATH=. python scripts/plot_fgc.py  ->  plots/fgc.png, ~/data/gc_fraction/fgc.csv
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import gcf

CODE = gcf.Path(__file__).resolve().parents[2]
mw_f = pd.read_csv(gcf.Path.home() / "data/mw_mdf/wmw/mass_by_mh_k1elfc.csv")
lmc_f = pd.read_csv(CODE / "lmc-ancient-disc/lmc-old-population/results/giant_present_mass.csv")
lmc_f = lmc_f[lmc_f.living_mass > 0]

mw = gcf.mw_gcs()
mw = mw[mw["IN_SITU"]]
omega = np.asarray(mw["NAME"]) == "NGC 5139"
lmc = pd.read_csv(gcf.DATA / "raw/lmc_fornax_gc_feh.csv").query("galaxy == 'LMC'")
mass = pd.read_csv(gcf.DATA / "raw/baumgardt_lg_masses.csv")
mass["name"] = mass.name.str.split().str.join(" ")
old = mass[(mass.galaxy == "LMC") & (mass.log_age > 9.9)]
lmc = old[["name", "mass"]].merge(lmc, on="name", how="left")
assert len(lmc) == 15 and lmc.feh.notna().all(), lmc[lmc.feh.isna()].name.tolist()
fnx_f = pd.read_csv(gcf.DATA / "fornax_mass.csv").query("living_mass > 0")
fnx = pd.read_csv(gcf.DATA / "raw/lmc_fornax_gc_feh.csv").query("galaxy == 'Fornax' and name != 'Fornax 6'")
fnx = fnx.merge(mass[mass.galaxy == "Fornax"][["name", "mass"]], on="name", how="left")
assert len(fnx) == 5 and fnx.mass.notna().all()


def binned(feh, m, edges):
    return np.histogram(feh, edges, weights=m)[0], np.histogram(feh, edges)[0]


mw_edges = np.r_[mw_f.mh_centre - 0.1, mw_f.mh_centre.iloc[-1] + 0.1]
mw_gc, mw_n = binned(mw["FEH"], mw["MASS"], mw_edges)
mw_gc_noom, _ = binned(mw["FEH"][~omega], mw["MASS"][~omega], mw_edges)
mini = 10 ** np.asarray(mw["MINI"])   # MINI = log10 Baumgardt initial mass
mw_gc_init, _ = binned(mw["FEH"], mini, mw_edges)
mw_gc_init_noom, _ = binned(mw["FEH"][~omega], mini[~omega], mw_edges)
lmc_edges = np.r_[lmc_f.feh_lo, lmc_f.feh_hi.iloc[-1]]
lmc_gc, lmc_n = binned(lmc.feh, lmc.mass, lmc_edges)
fnx_edges = np.r_[fnx_f.feh_lo, fnx_f.feh_hi.iloc[-1]]
fnx_gc, fnx_n = binned(fnx.feh, fnx.mass, fnx_edges)

rows = [dict(galaxy="MW", feh=x, n_gc=n, m_gc=gi, m_gc_no_omegacen=gi2, m_gc_present=g, m_living=f, f_living=gi / f,
             f_present_living=g / f)
        for x, n, g, gi, gi2, f in zip(mw_f.mh_centre, mw_n, mw_gc, mw_gc_init, mw_gc_init_noom, mw_f.M_total_Msun)]
rows += [dict(galaxy="LMC", feh=x, n_gc=n, m_gc=g, m_living=fl, f_living=g / fl, m_present=fp, f_present=g / fp,
              f_living_lo=g / r.living_mass_p84, f_living_hi=g / r.living_mass_p16)
         for (_, r), x, n, g, fl, fp in zip(lmc_f.iterrows(), lmc_f.feh, lmc_n, lmc_gc, lmc_f.living_mass, lmc_f.present_mass)]
rows += [dict(galaxy="Fornax", feh=x, n_gc=n, m_gc=g, m_living=r.living_mass, f_living=g / r.living_mass,
              f_living_lo=g / (r.formed_p84 * r.living_fraction), f_living_hi=g / (r.formed_p16 * r.living_fraction))
         for (_, r), x, n, g in zip(fnx_f.iterrows(), fnx_f.feh, fnx_n, fnx_gc)]
tab = pd.DataFrame(rows)
tab.to_csv(gcf.DATA / "fgc.csv", index=False)

fig, (a0, ax) = plt.subplots(1, 2, figsize=(13, 5))
for name, edges, field, gc, c in [("MW in-situ", mw_edges, mw_f.M_total_Msun, mw_gc_init, "tab:orange"),
                                  ("LMC", lmc_edges, lmc_f.living_mass, lmc_gc, "tab:blue"),
                                  ("Fornax", fnx_edges, fnx_f.living_mass, fnx_gc, "tab:green")]:
    e = np.array(edges, float)
    e[0] = max(e[0], e[1] - 0.5)          # open bin < -2 drawn as -2.5..-2.0
    e[-1] = min(e[-1], e[-2] + 0.5)
    field = np.asarray(field)
    mp = field[np.asarray(edges)[1:] <= -1.3 + 1e-9].sum()
    a0.stairs(field, e, color=c, lw=2, label=f"{name.split()[0]} field (living); M([Fe/H]<-1.3) = {mp:.1e}")
    a0.stairs(gc, e, color=c, lw=1.5, ls="--", label=f"{name} GCs" + (" (initial mass)" if name.startswith("MW") else ""))
a0.stairs(mw_gc, mw_edges, color="tab:orange", alpha=.5, lw=1, ls=":", label="MW in-situ GCs (present mass)")
a0.set(yscale="log", xlabel="[Fe/H]", ylabel=r"mass per bin [M$_\odot$]", xlim=(-2.5, 0.5), ylim=(1e4, 2e10))
a0.set_title("Mass per bin (MW 0.2 dex, LMC/Fornax 0.25 dex): field solid, GCs dashed", fontsize=10)
a0.legend(fontsize=7.5)
m, l, fx = (tab[(tab.galaxy == g) & (tab.n_gc > 0)] for g in ("MW", "LMC", "Fornax"))
ax.plot(m.feh, m.f_living, "o-", c="tab:orange", label="MW in-situ (BK): GC initial / living field")
ax.plot(m.feh, m.m_gc_no_omegacen / m.m_living, "o--", mfc="none", c="tab:orange", label=r"MW, initial, without $\omega$ Cen")
ax.plot(m.feh, m.f_present_living, ".:", c="tab:orange", alpha=.5, label="MW: GC present / living field")
ax.errorbar(l.feh, l.f_living, [l.f_living - l.f_living_lo, l.f_living_hi - l.f_living], fmt="s-", c="tab:blue",
            capsize=3, label="LMC, living field mass (p16–p84)")
ax.plot(l.feh, l.f_present, "s--", mfc="none", c="tab:blue", label="LMC, present field mass")
ax.errorbar(fx.feh, fx.f_living, [fx.f_living - fx.f_living_lo, fx.f_living_hi - fx.f_living], fmt="D", ms=7, c="tab:green",
            capsize=3, label="Fornax, living field mass (bootstrap p16–p84)")
for d in (m, l, fx):
    for x, y, n in zip(d.feh, d.f_living, d.n_gc):
        ax.annotate(str(n), (x, y), xytext=(4, 4), textcoords="offset points", fontsize=8)
ax.set(yscale="log", xlabel="[Fe/H]", ylabel=r"$f_{\rm GC} = M_{\rm GC}/M_*$", xlim=(-2.5, 0.5))
ax.set_title("GC mass fraction (MW: initial GC mass; LMC, Fornax: present; numbers = GCs per bin)", fontsize=9)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("plots/fgc.png", dpi=120)
print(tab.to_string(float_format="%.3g"))
print("MW in-situ GCs outside field bins:", int(((mw["FEH"] < mw_edges[0]) | (mw["FEH"] >= mw_edges[-1])).sum()))
