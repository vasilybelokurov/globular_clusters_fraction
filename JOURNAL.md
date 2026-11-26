# JOURNAL — GC mass fraction vs metallicity (globular_clusters_fraction)

Running log of what was done, decided and found. Newest entries at the bottom of each day.
Rule: updated by the agent at every milestone, result, decision, bug or agent report (see AGENTS.md).

---

## 2026-11-26 (Thu) — set-up

- Created `AGENTS.md` (least code that works; failing test first; plot before claiming) and this journal.
- Mission: f_GC([Fe/H]) = M_GC / M_* per [Fe/H] bin for MW, LMC, Fornax (maybe Sgr, GSE), with M_* from each galaxy's MDF in M☉.
- Inputs located (not yet read into code):
  - **MW:** `~/data/mw_mdf/wmw/mass_by_mh_k1elfc.csv` — 14 bins of [M/H], centres −2.2…+0.4, width 0.2 dex; **living** mass;
    total 3.436e10 M☉. Metal-poor bins small: −2.2: 3.1e7, −2.0: 6.0e7, −1.8: 8.7e7, −1.6: 1.67e8 M☉. 29–62 % of mass at
    [M/H] < −0.8 lies in the XP-blind wedge (extrapolated).
  - **LMC:** `../lmc-ancient-disc/lmc-old-population/results/giant_present_mass.csv` — [Fe/H] bins with edges
    −inf, −2.0, −1.75, −1.5, −1.25, −1.0, −0.75, −0.5, −0.25, 0, +inf; columns formed/living/present mass with p16/p84 and model range.
    Totals (`giant_present_mass.json`): formed 2.45e9, living 1.16e9, present 1.49e9 M☉; [Fe/H] < −1.5: formed 5.2e7 M☉.
- Open issues before any computation:
  1. Mass definition mismatch: MW table is living mass, LMC has formed/living/present. GC masses are present-day
     (or initial if corrected). Choose one definition and convert.
  2. Metallicity scale: MW in [M/H] (Andrae+23 XP), LMC in [Fe/H]; GC catalogues in [Fe/H]. Need [α/Fe] conversion or a common scale.
  3. Bins differ (MW 0.2 dex centred on even tenths; LMC 0.25 dex). Rebin onto a common grid.
  4. GC catalogues with individual masses + [Fe/H] for MW, LMC, Fornax, Sgr, GSE still to be chosen and verified.
  5. Fornax stellar-mass MDF: needs the LMC-style completeness correction — method to be ported from lmc-old-population.
  6. Sgr and GSE: total stellar mass and MDF in M☉ are less certain; decide later.

## 2026-11-26 (Thu) — GC catalogues: MW in-situ (BK), LMC/Fornax masses

- **MW in-situ = BK boundary** found in `~/Dropbox (Personal)/idl/work/globular_clusters/gc_catalog_create.pro`
  (E–Lz mask; potential MN disc 6.8e10 + Hernquist bulge 5e9 + NFW 1e12, c=19.5; Sun R0=8, (U,V,W)⊙=(−9.3, 251.5, 8.59)).
  Its output `FLAG` (1 = in-situ) is in `~/IoA Dropbox/.../data/catalogues/gc_catalog_updated.fits` (165 GCs; MASS = Baumgardt present mass, FEH, LZ, EN, PROG_MA).
- Ported the boundary to Python: `gcf.py::bk_energy_boundary`, `bk_in_situ`. Test `tests/test_gcf.py` (2 passed): recomputed flag = IDL FLAG for 165/165; boundary continuous.
- Plot `plots/mw_insitu.png` (`scripts/plot_mw_insitu.py`): **107 in-situ, ΣM = 2.57e7 M☉; 58 accreted, ΣM = 1.16e7 M☉** (present masses).
- ⚠ ωCen (NGC 5139, 3.94e6 M☉, [Fe/H] −1.53) is BK in-situ — 15 % of in-situ GC mass, in a bin where the MW field mass is ~2e8 M☉. Decision needed: keep or exclude (nucleus of an accreted dwarf, PROG_MA = G-E/Seq).
- **LMC/Fornax masses:** Baumgardt Local Group N-body table (https://people.smp.uq.edu.au/HolgerBaumgardt/globular/lgclusters/parameter.html),
  parsed → `~/data/gc_fraction/raw/baumgardt_lg_masses.csv` (411 clusters: LMC 330, SMC 74, Fornax 5, M31 2). Same mass scale as the MW masses.
  Old LMC (log age > 9.9): 15 GCs — NGC 1466, 1754, 1786, 1835, 1841, 1898, 1916, 1928, 1939, 2005, 2019, 2210, 2257, Hodge 11, Reticulum.
  Fornax 1–5 present; Fornax 6 absent.
- Rejected: `lmc-old-population/results/rr_mass_clusters.csv` (Cuevas-Otahola+2025; 10 clusters incl. SMC NGC 121/361; masses on a 0.1-dex grid).
- [Fe/H] for LMC/Fornax GCs: literature compilation running (agent) → `~/data/gc_fraction/raw/lmc_fornax_gc_feh.csv`.

## 2026-11-26 (Thu) — first f_GC([Fe/H]) for MW and LMC

- LMC/Fornax [Fe/H] compiled (agent) → `~/data/gc_fraction/raw/lmc_fornax_gc_feh.csv` (PDFs in `raw/literature/`).
  Preferred: Mucciarelli+2021 NatAs (homogeneous high-res, 11/15 LMC GCs; https://arxiv.org/abs/2110.10561); CaT for Reticulum, NGC 1841
  (Grocholski+2006), NGC 1928, 1939 (Piatti+2018, poorly determined). Fornax: Letarte+2006 (F1–3), Larsen+2018 (F4–5); Fornax 6 ~2 Gyr, not old.
  Spot-checked M21 Table 1 (NGC 1466, 1898, 2257): match. Johnson+2006 NGC 2005 −1.47 inconsistent with its own per-star values — not used.
- Bug: ωCen exclusion silently failed (FITS NAME bytes vs str under np.asarray). Failing test `test_mw_names_matchable` → fix
  `convert_bytestring_to_unicode()` in `gcf.mw_gcs` → 3 tests pass.
- `scripts/plot_fgc_mw_lmc.py` → `plots/fgc_mw_lmc.png`, table `~/data/gc_fraction/fgc_mw_lmc.csv`.
  f = M_GC(present, Baumgardt) / M_*(living). MW [M/H] treated as [Fe/H]; GCs on each galaxy's own field bins. LMC: 15 old GCs (log age > 9.9).
  - MW in-situ: 0.02–0.036 at [Fe/H] ≤ −1.6 (ωCen bin −1.6: 0.036 → 0.012 without), 0.006 at −1.4…−1.2, falling to 3e-4 at −0.2; no GCs at ≥ 0. 4 in-situ GCs below −2.3 outside field bins (dropped).
  - LMC: 0.11 (−2.25, 2 GCs), 0.025 (−1.88, 1), **0.21 (−1.62, 8 GCs)**, 0.018 (−1.38, 3), 0.003 (−1.12, 1); none at > −1.0. Present-mass denominator ×0.7.
  - At [Fe/H] < −1.5 LMC f_GC ≈ 3–7× MW; at −1.4…−1.1 comparable (0.003–0.018 vs 0.006).
- Caveats: GC present mass vs field living mass (MW biased high ×1.2–1.4); MW field error per bin not propagated; small-number noise (1–4 GCs/bin).

## 2026-11-26 (Thu) — Fornax added

- **Fornax field mass per [Fe/H]** `scripts/build_fornax_mass.py` → `~/data/gc_fraction/fornax_mass.csv`, `plots/fornax_mass.png`.
  LMC-style: Gaia DR3 parent (PM χ² < 9.21, parallax ≈ 0, PARSEC RGB box, 16.5 ≤ G0 ≤ 19.5; 5809 stars, background 22 from r_ell > 90′),
  Battaglia+2006 CaT members (VizieR J/A+A/459/423 via CfA mirror — CDS unreachable; 498/562 matched) weighted to the parent in r_ell × G0 cells;
  29 % of parent in cells without spectra (faint, r_ell < 20′) filled by rescaling within each r_ell bin. η from lmcold PARSEC grid + age windows A (ref), B, C.
  - Total **formed 3.42e7 M☉ (priors 2.76–4.12e7), living 1.54e7**. Test vs de Boer+2012 (2012A&A...544A..73D): 4.3e7 formed within r_ell < 0.8° → ours 20 % low (prior B 4 %). Pass (< 25 %).
  - Metal-poor field: [Fe/H] < −2: formed 2.45e6, living 9.9e5 (33 spec stars).
  - [Fe/H] scale check: K10 − B06 = −0.23 dex (σ 0.21) for 19 common stars, all metal-rich; metal-poor end untested. CaT may misplace [Fe/H] < −2 stars (Starkenburg+2010) — systematic.
- **f_GC with Fornax** `scripts/plot_fgc.py` (renamed from plot_fgc_mw_lmc.py) → `plots/fgc.png`, `~/data/gc_fraction/fgc.csv`.
  Fornax 1–5 (Baumgardt masses; [Fe/H] Letarte+2006 / Larsen+2018); Fornax 6 excluded (~2 Gyr). Baumgardt names have double spaces — normalised.
  - Fornax [Fe/H] < −2 (F1, 2, 3, 5; 2.93e5 M☉): **f = 0.30 (0.23–0.43)** of living field mass (0.12 of formed). [Fe/H] −1.5…−1.25 (F4): 0.07 (0.05–0.11).
  - Ordering at [Fe/H] < −2: Fornax 0.30 > LMC 0.11 > MW 0.01–0.02.

## 2026-11-26 (Thu) — f_GC figure made 2-panel
- `plots/fgc.png` now: (left) dM/d[Fe/H] [M☉/dex] of field (living, solid) and GCs (dashed) for MW (GCs in-situ only), LMC, Fornax;
  open bin < −2 drawn as −2.5…−2.0; (right) f_GC as before. Per-dex normalisation because bin widths differ (MW 0.2, LMC/Fornax 0.25).

## 2026-11-26 (Thu) — MW field MDF vs Kurbatov Aurora (user query)
- User: Kurbatov reports ~0.5e9 M☉ for Aurora, the plot showed much more. Cause: left panel was dM/d[Fe/H] (per dex); MW bins are 0.2 dex → ×5 visual.
  Table itself: M([M/H] < −1.3) = 5.53e8, < −1.1: 9.3e8, < −0.9: 1.29e9 (k1elfc). mw_mdf earlier GMRF check: 4.5e8 vs paper 4.3e8 (−2.9…−1.3); k1elfc is 28 % above 4.3e8.
- Left panel switched to mass per bin, legend gives M([Fe/H] < −1.3) per galaxy (MW 5.5e8, LMC 2.1e7, Fornax 2.8e6).
- **Open issue (real):** MW field table = all stars (in-situ + accreted halo), GC numerator = in-situ only → MW f_GC biased low where the accreted halo matters ([Fe/H] ≈ −1.5…−1).
  The Kurbatov pipeline (`../mw_mdf/gaiaunlimited-aurora/Aurora`, runs pl-sc-* with a GSE component vs pl-nogse-*) may give an Aurora-only (in-situ) mass per [M/H] bin — to be checked.

## 2026-11-26 (Thu) — MW curve with GC initial masses
- `MINI` in gc_catalog_updated.fits = log10 Baumgardt initial mass (identical to `../globular_clusters_imf` initial_mass_msun for 149/165 matched; max |Δ| = 0).
- In-situ GCs: Σ present 2.57e7, Σ initial 1.65e8 M☉ (×6.4). Added to `plots/fgc.png` (red dotted, both panels; column f_initial_living in fgc.csv).
- f_initial/living: ≈ 0.07–0.11 at [Fe/H] −2.2…−1.4 (flat; ωCen bin no longer stands out), 0.045–0.05 at −1.2…−1.0, ~0.007–0.01 at −0.8…−0.4, 1.5e-3 at −0.2.
- Note: initial GC mass vs living field mass mixes definitions; formed field mass (×~2.4 at low [Fe/H], LMC living fractions 0.40–0.43) would lower it to ≈ 0.03–0.045.

## 2026-11-26 (Thu) — MW main curve = GC initial mass (user decision)
- USER DECISION: MW f_GC uses in-situ GC **initial** masses (Baumgardt MINI). `plots/fgc.png`: solid orange = initial/living field; dashed = initial without ωCen;
  faint dotted = present mass (reference). fgc.csv MW columns: m_gc = initial, m_gc_present, f_living (initial-based), f_present_living.
- Without ωCen (initial): bin −1.6 0.108 → 0.061.
- Inconsistency now: LMC/Fornax still present GC masses (Baumgardt LG table gives no initial masses) → MW vs dwarfs not like-for-like.

## 2026-11-26 (Thu) — STATUS at session end (resume here)

**Reproduce** (venv `~/Work/venvs/.venv`, from project root):
`PYTHONPATH=. python -m pytest -q tests` (3 pass) · `PYTHONPATH=. python scripts/plot_mw_insitu.py` · `PYTHONPATH=. python scripts/build_fornax_mass.py` · `PYTHONPATH=. python scripts/plot_fgc.py`

**Files:** code `gcf.py` (BK boundary, MW GC loader), `scripts/`, `tests/`; plots `plots/{mw_insitu,fornax_mass,fgc}.png`.
Data `~/data/gc_fraction/`: `fgc.csv` (main result table), `fornax_mass.csv`, `fornax_gaia.parquet` (cached WSDB query);
`raw/`: `baumgardt_lg_masses.csv` (+ html source), `lmc_fornax_gc_feh.csv`, `literature/` (23 PDFs), `fornax/` (B06, K10 VizieR tables).
Removed superseded `fgc_mw_lmc.csv`. Project is **not** under git yet.

**Current headline (plots/fgc.png):** MW in-situ (BK) GC initial / living field: 0.07–0.11 at [Fe/H] ≤ −1.4, 0.045–0.05 at −1.2…−1.0, ~0.007–0.01 at −0.8…−0.4.
LMC (present GC / living field): 0.11, 0.025, 0.21, 0.018, 0.003 for bins −2.25…−1.12. Fornax: 0.30 (< −2), 0.07 (−1.38).

**Open decisions / next steps:**
1. ωCen in or out of MW in-situ (BK says in-situ).
2. Mass definitions: MW uses GC initial mass, LMC/Fornax GC present; field = living for all. Options: formed field mass (LMC/Fornax tables have it; MW needs an AMR-based living→formed conversion); initial masses for LMC/Fornax GCs (not in Baumgardt LG table).
3. MW denominator includes accreted halo; check whether the Kurbatov Aurora pipeline (`../mw_mdf/gaiaunlimited-aurora/Aurora`, pl-sc-* with GSE component) gives in-situ mass per [M/H] bin.
4. Common [Fe/H] grid across galaxies; [M/H]≈[Fe/H] assumption for MW field; Fornax CaT scale at [Fe/H] < −2 untested.
5. Sgr and GSE (PROG_MA column in gc_catalog_updated.fits gives Massari+2019 progenitor labels).
6. Put the project under git.
