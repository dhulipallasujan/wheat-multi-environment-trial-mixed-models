# California Wheat Trial Analysis

## In plain language

This project analyzes real California wheat field-trial records from 2005–2018: 13,953 usable plot yields for 211 entries across 103 year-location environments. In the multi-environment model, Entry 1845 had the highest estimated mean yield (7,784 kg/ha), but that ranking is uncertain and the entries were not all tested in the same number of environments. The estimated entry-by-environment variation was substantial (SD about 717 kg/ha), so an entry's performance can change across environments. In the separate COLUSA 2010 trial, Entry 1478 had the highest observed mean (10,306 kg/ha), but Tukey comparisons did not distinguish it from several other top entries.

## Research questions

1. Did wheat entries differ in yield in the balanced COLUSA 2010 trial after accounting for blocks?
2. Which entries had higher estimated mean yield across the available year-location environments?
3. How much did entry rankings vary among environments?

## Data and provenance

The plot-level data include entry number, year, location, block, and yield.

- Dataset and source: [Dryad dataset](https://doi.org/10.5061/dryad.bf8rt6b)
- Study: George and Lundy (2019), [Crop Science paper](https://doi.org/10.2135/cropsci2019.01.0010)
- The Dryad usage notes state that yield is measured in kg/ha and that the dataset is released under CC0. The CSV itself has no separate units column.
- The downloaded file is `data/George_et_al_2019.csv`. The 43 `#VALUE!` yield cells are treated as missing; they are not replaced or imputed.

## Analysis

Open `wheat-trial-project.Rproj` in RStudio, then run `source("analysis.R")` from the project root. The script creates the figures and result files listed below. Required packages are `lme4`, `emmeans`, and `multcompView`; install them once in R with:

```r
install.packages(c("lme4", "emmeans", "multcompView"))
```

The script cleans the yield column, analyzes COLUSA 2010, and fits a multi-environment mixed model. In the multi-environment model, entry is treated as a fixed effect, while year-location environment, block within environment, and entry-by-environment deviations are random effects. Model-adjusted entry means come from `emmeans`. The ranking summary also reports how many environments each entry appeared in and the standard deviation of its observed within-environment ranks. These are descriptive rank summaries, not a formal stability test such as AMMI or Finlay-Wilkinson regression.

The COLUSA 2010 analysis uses an ANOVA with entry and block terms, Tukey-adjusted comparisons, model-based 95% confidence intervals, and diagnostic flags. Its chart shows every one of the 29 entries. Points are means, horizontal lines are 95% confidence intervals, and shared Tukey letters indicate pairs that were not significantly different at the 5% level. Shared letters do not prove entries are identical.

## Results and interpretation

### Multi-environment analysis

The cleaned data contain 13,953 plot records, 211 entries, and 103 year-location environments. Entry 1845 ranked first by model-adjusted mean (7,784 kg/ha; asymptotic 95% CI: 7,112–8,456), followed by Entry 1847 (7,505 kg/ha; tested in 13 environments) and Entry 1843 (7,453 kg/ha; tested in 7 environments). Entry 1845 appeared in 7 environments. The estimated standard deviation for entry-by-environment deviations was about 717 kg/ha, indicating meaningful variation in performance across environments.

This is an unbalanced dataset: entries were tested in different numbers and combinations of environments. The model-based ranking is a useful summary of these data, not proof that the first-ranked entry is universally best. The estimated environment standard deviation was about 1,637 kg/ha, larger than the entry-by-environment SD, reflecting strong differences among trial environments.

### COLUSA 2010 example

This subset had 29 entries and 29 records in each of four blocks. The overall entry effect was significant (F = 16.419, p < 2 × 10⁻¹⁶), as was the block effect (p = 0.0284). Entry 1478 had the highest observed mean yield, 10,306 kg/ha. Tukey-adjusted comparisons did not distinguish it from the next four entries by mean (adjusted p-values > 0.94). The overall ANOVA establishes that at least one entry mean differs; it does not show that Entry 1478 is better than every other entry.

## Limitations and next steps

The mixed model summarizes an unbalanced set of trials and uses a simple random-intercept structure. Its rankings should be interpreted with the environment coverage and uncertainty in mind. The dataset provided here has no genetic-background information, so these results cannot explain why entries respond differently. The current rank summaries are exploratory; formal stability analysis and a sensitivity analysis of flagged records would be useful next steps. Diagnostic plots are saved for review; flagged observations should be checked against the source before any decision to exclude them.

## Project files

- `analysis.R` — cleaning, single-environment and multi-environment analyses
- `LICENSE` — MIT license for the analysis code
- `data/George_et_al_2019.csv` — downloaded source data
- `figures/colusa_2010_all_entries.png` — all 29 COLUSA 2010 entries with confidence intervals and Tukey groups
- `figures/model_diagnostics.pdf` — COLUSA 2010 model diagnostic plots
- `figures/multi_environment_diagnostics.pdf` — multi-environment model residual plot
- `results/colusa_2010_entry_means.csv` — COLUSA 2010 means, intervals, and Tukey groups
- `results/tukey_colusa_2010.csv` — Tukey-adjusted pairwise comparisons
- `results/multi_environment_entry_means.csv` — model-adjusted means, intervals, coverage, and rank variation
- `results/entry_environment_means.csv` — observed entry means and ranks by environment
- `results/multi_environment_variance_components.csv` — estimated random-effect variances
- `results/observations_to_review.csv` — records flagged for diagnostic review
- `results/session_info.txt` — R and package versions used for the run
