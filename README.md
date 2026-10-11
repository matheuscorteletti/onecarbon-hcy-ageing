# One-carbon metabolism across age strata — R pipeline (Clinical Nutrition)

Figures and tables for the manuscript testing **H1**: in hospitalised adults living under mandatory
folic-acid fortification, the hyperhomocysteinaemia of the oldest old (≥85 y) is explained mainly by
kidney function and active B12 (not folate), with an independent contribution of malnutrition.

## Run

**Option A — single file (recommended).** Open `onecarbon_all_in_one.R` in RStudio and click **Source**.
A window asks for the Excel workbook; the data are read once and every figure and table is written to
`onecarbon_output/` next to the workbook (`figures/`, `supplementary/`, `tables/`). Missing packages are
installed automatically. Edit `N_BOOT_USER` (1000 final / 200 draft) and `SAVE_PANELS_USER` at the top.

**Option B — modular scripts** (same code split by figure; useful for editing one figure):

```r
source("install_packages.R")   # once
# put the workbook in the project root (or set ONECARBON_DATA)
Rscript run_all.R              # ~6 min with 1000 bootstrap replicates
```

| Variable | Default | Meaning |
|---|---|---|
| `ONECARBON_DATA` | `banco_de_dados_onecarbon.xlsx` | path to the workbook (sheets *Data*, *Data_Dictionary*, *Genotype_QC*) |
| `ONECARBON_NBOOT` | `1000` | bootstrap replicates (use `200` for quick drafts) |
| `ONECARBON_PANELS` | `FALSE` | `TRUE` also exports every panel as its own TIFF (for re-layout in Illustrator/BioRender) |

Each script is independent (`source("R/fig4_determinants_hcy.R")` works on its own).

## Output

* `output/figures/Figure1…7_*.tiff` — 180 mm wide (double column), 600 dpi, LZW, RGB. **No panel letters
  and no captions** (to be added at layout). Font: Arial (Liberation Sans, metric-identical, if Arial is absent).
* `output/supplementary/FigS1…S5_*.tiff`
* `output/tables/Tables_1-3.docx`, `Supplementary_Tables.docx` (+ CSV of every table and of the numbers
  quoted in the text: Oaxaca, LMG, g-computation, mediation, SEM paths, CV-AUC, variant associations).

Visual system: one fixed colour per group in every figure (Community 56–75 grey · Hospital 56–75 blue ·
76–84 amber · ≥85 crimson), colour-blind safe (distinct in hue and lightness); diverging blue–white–red
for β/ρ; sequential cream→dark red for Hcy; 7 pt base text.

## Script structure

| Script | Content |
|---|---|
| `R/00_setup.R` | packages, data preparation, derived variables, palette, theme, export |
| `R/01_models.R` | complete-case H1 set, hierarchical blocks M1–M5, g-computation engine, rank-preserving intervention |
| `R/fig1…fig7_*.R` | one script per figure |
| `R/tables_main.R` | Tables 1–3 (flextable → Word) |
| `R/supplementary.R` | Figs S1–S5, Tables S1–S4 |

## Panel map (in reading order inside each composite)

**Figure 1 — Design and causal framework (5 panels).** Flow 4 groups → analytic sets · DAG (bold = H1 paths) ·
age ridgelines with strata cut-points · data-availability matrix · % of life under fortification.

**Figure 2 — Nutritional phenotype, Q1 (11 panels).** GLIM / MNA-SF / EWGSOP2 prevalence · alluvial
age→MNA-SF→GLIM→EWGSOP2 · phase angle raincloud · handgrip by sex with EWGSOP2 cut-offs · phase angle × ASMI
(GLIM colour, 68% ellipses per stratum) · energy × protein intake vs ESPEN targets · PCA loadings · PCA scores
with centroid trajectory · Gaussian-mixture phenotypes (profile heatmap; data-driven names Robust / Acute
weight loss / Sarcopenic / Cachectic) · phenotype composition by group · Hcy by phenotype within strata.

**Figure 3 — Diet to status under fortification, Q2 (10 panels).** Intake/EAR rainclouds (5 nutrients) ·
mirror plot intake < EAR vs biochemical insufficiency · folate sources (natural / fortified flour / supplements,
% > UL) · folic acid → UMFA by DHFR genotype · 4 intake→biomarker scatters (ρ overall and per group) ·
diet–biomarker ρ with bootstrap CI per group · RBC folate × holoTC quadrants coloured by Hcy.

**Figure 4 — Determinants of Hcy, Q3 (12 panels; core of the paper).** Hcy by group (% > 15) · age gradient of
determinants (% of 56–75 y) · attenuation M1→M5 · Oaxaca–Blinder decomposition (pooled, bootstrap CI) ·
LMG relative importance (bootstrap CI) · adjusted GAM curves by stratum for eGFR-cys, holoTC and RBC folate
(dashed = crude folate curve) · **folate paradox** (crude vs partial ρ) · joint eGFR × holoTC surface with
stratum centroids · MMA vs eGFR by holoTC (why holoTC) · quantile regression of the ≥85 y effect, M1 vs M5.

**Figure 5 — Counterfactual, Q4 (6 panels).** g-computation for 76–84 y and ≥85 y under rank-preserving
interventions (younger folate [negative control], eGFR, holoTC, both, + nutrition) · % Hcy > 15 under each
scenario · dose–response of the share of the gap closed · interventional mediation waterfall (kidney, B12,
nutrition, overlap, direct) · predicted Hcy distributions · robustness across 9 sensitivity sets.

**Figure 6 — Genetics vs modifiable determinants, Q5 (7 panels).** LMG variance partition by group · per-SD
effects (GRS vs eGFR, holoTC, folate, phase angle) by group · forest of the 24 variants (FDR; SHMT1 null control)
· variant × biomarker heatmap (functional validation: TCN2→holoTC, FUT2→B12, DHFR→UMFA…) · MTHFR 677 × RBC folate ·
MTHFR 677 × EGRAC · GRS tertile × low holoTC.

**Figure 7 — Integration and translation, Q6 (7 panels).** EBIC-glasso partial-correlation network (edge alpha =
bootstrap inclusion) · clustered partial-Spearman heatmap (one-carbon metabolites × nutritional markers) · SEM
(malnutrition → B12 → Hcy → SAM/SAH → 5-mC; Hcy → ADMA/GDF-15; fit indices) · Hcy vs SAM/SAH, 5-mC, ADMA, GDF-15 ·
conditional inference tree for Hcy > 15 · etiological profile of hyperhomocysteinaemia by stratum · cross-validated
AUC of nested determinant sets.

**Supplementary.** S1 folate "ceiling" (H1b; piecewise slopes) · S2 folate × B12 (H1c) for Hcy and MMA ·
S3 all 144 gene × nutrient interactions · S4 genotype QC (HWE, call rate, EAF) · S5 exploratory mortality.
Tables S1 sensitivity of M5 (incl. multiple imputation, eGFR-creatinine, original design without 76–84 y) ·
S2 Tobit vs OLS for serum folate (ULOQ 24 ng/mL) · S3 genotype QC · S4 all GxN interactions.

## Notes for the Methods section

* Primary set: hospital participants with valid Hcy (3 strata); complete cases for models; MI (m = 20) as sensitivity.
* Outcome log-transformed; effects reported as % change = 100 × (e^β − 1); back-transformation with Duan smearing.
* Interventions = rank-preserving distributional shift: each older participant receives the 56–75 y quantile
  matching their own rank (interventional effects; bootstrap stratified by age group).
* Levodopa users (n = 5 in the model set) are handled by exclusion (sensitivity), not adjustment.
* Packages that the original plan listed but are replaced here by transparent code: `oaxaca` (pooled
  twofold decomposition coded directly — identical point estimates), `CMAverse` (interventional mediation via
  g-computation), `dagitty/ggdag` (DAG drawn in ggplot2), `forestploter` (ggplot2), `gtsummary` (flextable),
  `ComplexHeatmap` (ggplot2 tiles with hclust order, so panels compose with patchwork).

## Repository and data policy

* **The data workbook is never committed** (`*.xlsx`, `*.csv` and all generated outputs are listed in
  `.gitignore`). Participant-level data stay on the institutional drive under the ethics approval.
* Keep the repository **private** until the paper is accepted; make it public (and add a licence, e.g. MIT)
  at publication and cite the release DOI (Zenodo–GitHub integration) in the *Data and code availability* statement.

### Putting the code on GitHub

*Web upload (no command line)*
1. github.com → **New repository** → name `onecarbon-hcy-ageing`, visibility **Private**, README **Off**,
   .gitignore **No .gitignore**, licence **No license** → **Create repository**.
2. On the empty repository page click **uploading an existing file**.
3. Drag in `README.md`, `.gitignore`, `onecarbon_all_in_one.R`, `run_all.R`, `install_packages.R` and the
   folder `R/` (on macOS press Cmd+Shift+. in Finder to see `.gitignore`). Do **not** drag the `.xlsx` or `output/`.
4. Commit message `Initial pipeline` → **Commit changes**.

*Command line (Terminal, inside the project folder)*
```bash
git init
git add .                      # .gitignore keeps data and outputs out
git status                     # check: no .xlsx, no output/
git commit -m "Initial pipeline"
git branch -M main
git remote add origin https://github.com/matheuscorteletti/onecarbon-hcy-ageing.git
git push -u origin main
```

*Later changes*: `git add -A && git commit -m "describe the change" && git push`
(or edit/upload files on the website).
