# n-alkane vegetation correction for δ²H of precipitation

An R script that reconstructs the hydrogen-isotope composition of precipitation
(δ²H<sub>prc</sub>) from sedimentary *n*-C<sub>29</sub> alkane δ²H. It corrects for
changes in the grass vs. woody-plant contribution and propagates the uncertainty
with a Monte Carlo simulation, following **Santos et al. (2026)**.

It also fits GAM smooths, finds where the change is statistically significant,
and writes Grapher/Excel-ready tables, publication figures and a methods summary.

## Method

| Step | Equation |
|---|---|
| Grass fraction (relative-abundance index; Schäfer et al., 2016) | f<sub>GR</sub> = (C31 + C33) / (C27 + C31 + C33) |
| Mixed apparent fractionation | ε<sub>mix</sub> = f<sub>GR</sub>·ε<sub>GR</sub> + (1 − f<sub>GR</sub>)·ε<sub>WP</sub> |
| Precipitation δ²H | δ²H<sub>prc</sub> = (δ²H<sub>C29</sub> + 1000) / (1 + ε<sub>mix</sub>/1000) − 1000 |

* ε<sub>WP</sub> = −110 ± 21 ‰ and ε<sub>GR</sub> = −165 ± 25 ‰ (normal distributions). You can change both in the settings.
* 10,000 Monte Carlo draws per sample. They include endmember uncertainty and, when you give it, the analytical SD of δ²H<sub>C29</sub>.
* **n-C33 is optional.** With `missing_C33_policy = "zero"` (the default), a missing C33 value, or a missing C33 column, gives f<sub>GR</sub> = C31 / (C27 + C31). With `"exclude"`, those samples are dropped. Figures and outputs always show which formula was used.
* GAMs use `mgcv::gam`: Gaussian family, thin-plate spline, REML. The basis size k is set from the data and checked with `k.check`. The script fits GAMs to raw δ²H<sub>C29</sub>, f<sub>GR</sub> and δ²H<sub>prc</sub>.
* **Significant change.** The script computes the first derivative of each GAM with a pointwise 95 % CI. Where the whole CI is above or below zero, the change is significant (Simpson, 2018). Rates are always given **forward in time**, so "increase" means the value gets higher towards the present.

## Quick start

1. Install R (≥ 4.0) and the packages:
   ```r
   install.packages(c("mgcv", "readxl", "writexl"))
   ```
2. Clone or download this repository and run the example:
   ```bash
   Rscript run_nalk_vegetation_correction.R
   ```
   or open the script in RStudio and click **Source**.
3. To use your own data, edit **section 1 (USER SETTINGS)** at the top of the script. Usually you only change `input_file` and, if needed, the column names and age settings. You can also pass the file path on the command line:
   ```bash
   Rscript run_nalk_vegetation_correction.R path/to/my_data.xlsx "n-alkanes"
   ```

## Input

An `.xlsx`, `.xls` or `.csv` table with one row per sample. For Excel files, the script reads the sheet named `n-alkanes` if there is one, otherwise the first sheet. The default column names are below; you can rename any of them in the settings.

| Column (default name) | Required | Description |
|---|---|---|
| `age_ka_bp` | ✔ | Sample age. cal ka BP, cal yr BP and year CE all work (set `age_is_BP`, `age_unit`, `age_axis_label`) |
| `nalk_C29_d2H_permille` | ✔ | *n*-C<sub>29</sub> δ²H (‰ VSMOW) |
| `nalk_C27_conc_ug_g_dw`, `nalk_C31_conc_ug_g_dw` | ✔ | Concentrations, in any consistent unit |
| `nalk_C33_conc_ug_g_dw` | – | Concentration. If missing, handled by `missing_C33_policy` |
| `nalk_C29_d2H_sd` | – | Analytical 1 SD (‰). If missing, treated as 0 and flagged |
| `nalk_C29_conc_ug_g_dw` | – | Carried through to the outputs only |
| `CPI_C23_33` | – | Samples with CPI below `cpi_review_threshold` (3) are flagged |
| `Sample_ID`, `depth_mid_cm` | – | Any `id_columns` are carried through to every output |
| `age_ka_bp_q025`, `_q16`, `_q84`, `_q975` | – | Age-model quantiles. They enable the age-uncertainty envelope |

`example/example_input.xlsx` (and `.csv`) contains **synthetic data** in the expected format. Use it as a template.

## Outputs

Each input file gets its own folder:

```
outputs/<input-file-name>/
├── <name>_RESULTS.xlsx          main workbook (open the README sheet first)
├── <name>_METHODS_SUMMARY.txt   methods text, significant intervals, diagnostics, caveats
├── 01_tables/                   CSV copies: observations, GAM curves, significant
│                                intervals, Monte Carlo summary, age envelope
├── 02_figures/png, pdf          3-panel figure, δ²Hprc, rates of change, age envelope
└── 03_diagnostics/              GAM diagnostics, residual plots, QC summary,
                                 excluded rows, run settings, run manifest, sessionInfo
```

**Workbook sheets:** README · Observations · GAM · Significant_intervals · MC_summary ·
GAM_diagnostics · QC_summary · Excluded_rows · Age_envelope (only if age quantiles are given) · Run_info

### Plotting significant parts of a curve (Grapher, Excel, …)

The **GAM** sheet has the following columns for each series (`d2H_C29`, `f_GR`, `d2Hprc`):

| Column | Content |
|---|---|
| `<s>_GAM`, `_lower_95`, `_upper_95` | Fitted curve and 95 % confidence band |
| `<s>_GAM_signif` | The fitted value **only where the change is significant**; blank elsewhere |
| `<s>_GAM_signif_increase` / `_decrease` | The same values, split by direction (forward in time) |
| `<s>_rate_per_<unit>` (+ CI) | Rate of change |
| `<s>_signif_code` | +1 increase, −1 decrease, 0 not significant |
| `Notes` | Plain-text list of which series change significantly at that age, with rate and CI |

Plot `<s>_GAM` as a thin line. Then plot `<s>_GAM_signif` (or the increase and decrease columns) on top as a thick line. The blank cells break the line, so only the significant segments show. The **Significant_intervals** sheet summarises every segment: start and end age, net change, mean and peak rate, and the number of samples inside. Its *caution* column marks segments with fewer than 3 samples or at the edge of the record.

## Limitations and appropriate use

**Relative-abundance (R.A.) index.** f<sub>GR</sub> assumes that *n*-C<sub>31</sub> and *n*-C<sub>33</sub> come mainly from grasses and *n*-C<sub>27</sub> mainly from woody plants. This does not hold for every species. As discussed in Santos et al. (2026), some trees, such as *Fraxinus* and *Acer*, also produce substantial *n*-C<sub>31</sub> and *n*-C<sub>33</sub>. Species with high *n*-alkane production can therefore bias f<sub>GR</sub>, even when they make up a small part of the vegetation. The grass vs. woody-plant estimates from *n*-alkanes describe the wax source. They are **not a calibrated estimate of grass cover in the catchment** and cannot be used in place of pollen-based vegetation reconstructions.

**Transferability.** Vegetation-correction approaches depend on their calibration datasets and may not fully remove vegetation effects. Their transferability between records should be tested before they are applied in other regions and vegetation settings. Future work should test similar corrections in other regions and settings.

**Constant ε<sub>app</sub> endmembers.** Using a constant ε<sub>app</sub> for each vegetation endmember is a transparent, first-order way to assess the direction of vegetation effects and to approximate their magnitude. However, ε<sub>app</sub> may not be constant. The Monte Carlo simulation, with conservative spreads around each ε<sub>app</sub>, accounts for part of this. Residual physiological variability remains an additional source of uncertainty.

**Statistics.** The GAM significance flags are pointwise and exploratory. They are not simultaneous intervals and do not include chronology uncertainty. Check the residual autocorrelation in `03_diagnostics` before interpreting short-lived changes.

## Citation

If you use this script, please cite:

> Santos, R. N., Nelson, D. B., Klatt, A., Schubert, C. J., Dubois, N., De Jonge, C., & Ladd, S. N. (2026). Central European hydroclimate since the Younger Dryas inferred from vegetation-corrected sedimentary plant wax δ²H values. *Paleoceanography and Paleoclimatology*, 41, e2025PA005401. https://doi.org/10.1029/2025PA005401

See also `CITATION.cff` in the repository root.

## References

* Schäfer, I. K., et al. (2016). *SOIL*, 2, 551–564. https://doi.org/10.5194/soil-2-551-2016
* Sachse, D., et al. (2012). *Annu. Rev. Earth Planet. Sci.*, 40, 221–249. https://doi.org/10.1146/annurev-earth-042711-105535
* Simpson, G. L. (2018). *Front. Ecol. Evol.*, 6, 149. https://doi.org/10.3389/fevo.2018.00149
* Wood, S. N. (2017). *Generalized Additive Models: An Introduction with R* (2nd ed.). https://doi.org/10.1201/9781315370279

## License

See the repository root.
