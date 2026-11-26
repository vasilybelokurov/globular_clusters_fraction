"""E-Lz with the BK boundary, and mass vs [Fe/H] for MW in-situ / accreted GCs.

Usage: PYTHONPATH=. python scripts/plot_mw_insitu.py  ->  plots/mw_insitu.png
"""
import matplotlib.pyplot as plt
import numpy as np

import gcf

t = gcf.mw_gcs()
ins = np.asarray(t["IN_SITU"])
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

lz = np.linspace(-4.4, 4.4, 500)
ax[0].plot(lz, gcf.bk_energy_boundary(lz), "k-", lw=1)
for m, c, lab in [(ins, "tab:orange", f"in-situ ({ins.sum()})"), (~ins, "tab:blue", f"accreted ({(~ins).sum()})")]:
    ax[0].scatter(t["LZ"][m] * 1e-3, t["EN"][m] * 1e-5, s=12, c=c, label=lab)
    ax[1].scatter(t["FEH"][m], t["MASS"][m], s=12, c=c, label=lab)
ax[0].set(xlabel=r"$L_z$ [$10^3$ kpc km/s]", ylabel=r"$E$ [$10^5$ km$^2$/s$^2$]", xlim=(-4.5, 4.5), ylim=(-2.6, 0.2))
ax[1].set(xlabel="[Fe/H]", ylabel=r"present mass [M$_\odot$]", yscale="log", xlim=(-2.6, 0.2))
ax[0].legend()
fig.suptitle(f"MW GCs, BK boundary; in-situ mass {t['MASS'][ins].sum():.2e}, accreted {t['MASS'][~ins].sum():.2e} M$_\\odot$")
fig.tight_layout()
fig.savefig("plots/mw_insitu.png", dpi=120)
print("in-situ", ins.sum(), "mass", t["MASS"][ins].sum(), "| accreted", (~ins).sum(), t["MASS"][~ins].sum())
print("in-situ with FEH missing/NaN:", int(np.sum(~np.isfinite(t["FEH"][ins]))))
