# AGENTS.md

## Rules

- Write least code that works
- Always reproduce bugs with a failing test
- Show me a plot before claiming a result
- Update JOURNAL.md automatically at every milestone, result, decision, bug and agent report (dated, concise, file paths + key numbers). Newest entries at the bottom.
- Activate the venv before running anything: `source ~/Work/venvs/.venv/bin/activate`
- Every result comes from a script (`scripts/`); plots go to `plots/` (PNG unless PDF requested)
- Large data live in `~/data/`, never in the project folder; read sibling-project products in place, do not copy them

## Mission

Measure the mass fraction of stars in globular clusters as a function of metallicity:

    f_GC([Fe/H]) = M_GC([Fe/H]) / M_*([Fe/H])

per galaxy, where M_GC is the summed mass of the galaxy's GCs in a [Fe/H] bin and M_* is the galaxy's
stellar mass in the same bin from its MDF in absolute solar masses.

Galaxies: MW, LMC, Fornax; possibly Sgr and GSE.

## Inputs (sibling projects)

- MW MDF in M☉: `../mw_mdf` (headline per-[M/H] table `~/data/mw_mdf/wmw/mass_by_mh_k1elfc.csv`, living stars)
- LMC MDF in M☉: `../lmc-ancient-disc/lmc-old-population/results/giant_present_mass.csv` (formed / living / present)
- Fornax field: `scripts/build_fornax_mass.py` (Gaia DR3 parent + Battaglia+2006 CaT, LMC-style eta from lmcold)
- LMC/Fornax GC masses: Baumgardt Local Group N-body table (`~/data/gc_fraction/raw/baumgardt_lg_masses.csv`); [Fe/H]: `raw/lmc_fornax_gc_feh.csv`
- MW GCs: `gc_catalog_updated.fits` (Baumgardt masses incl. initial MINI, [Fe/H], E, Lz); in-situ = BK E–Lz boundary (`gcf.bk_in_situ`)
- GC catalogues and related work: `../globular_clusters_mw`, `../globular_clusters_in_dwarfs`, `../globular_clusters_completeness`

## Consistency requirements

- Numerator and denominator in the same mass definition (formed, present, or living) — state which.
- Same metallicity scale and bin edges for GCs and field stars ([Fe/H] vs [M/H]; state the conversion).
- GC masses: present-day vs initial — state which.
