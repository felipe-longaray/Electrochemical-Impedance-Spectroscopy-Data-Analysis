# EIS Equivalent Circuit Statistical Pipeline

[![R](https://img.shields.io/badge/R-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains a robust statistical pipeline for analyzing Electrochemical Impedance Spectroscopy (EIS) data. Specifically, it processes the individual parameters extracted from Equivalent Electrical Circuit (EEC) fitting (e.g., $R_{ct}$, $C_{dl}$, $R_s$) across complex, multi-variable experimental designs.

## Scientific & Statistical Context
In interfacial electrochemistry and biosensor development, evaluating surface modifications requires modeling the impedance response. Key parameters include:
* **$R_{ct}$ (Charge Transfer Resistance):** Dictates the kinetics of the Faradaic process at the electrode-electrolyte interface.
* **$C_{dl}$ (Double-Layer Capacitance):** Tracks changes in the dielectric nature and thickness of the interfacial layer.

When analyzing modified surfaces across multiple conditions (e.g., distinct surface chemistries, washing protocols, and plasma treatments), simple pairwise t-tests are statistically invalid. This script solves this by implementing:
1. **Factorial (3-Way) ANOVA:** To assess the main effects and interaction terms of the experimental matrix.
2. **Outlier Diagnostics:** Programmatically calculates **Cook's Distance** for all data points to objectively identify high-leverage outliers prior to contrast testing.
3. **Orthogonal Contrasts:** Utilizes custom contrast matrices via `emmeans` to test highly specific scientific hypotheses (e.g., comparing specific functionalized states).
4. **Holm-Bonferroni Corrections:** Strictly adjusts $p$-values to control the Family-Wise Error Rate (FWER) during multiple comparisons.

## Code Architecture
The code is written following the DRY (Don't Repeat Yourself) principle. A core functional loop automatically ingests the data, dynamically generates the ANOVA formulas for all five circuit parameters ($C_a, C_{dl}, R_a, R_{ct}, R_s$), evaluates Cook's distances, and executes the adjusted contrasts.

## Dependencies
* `dplyr`
* `stats`
* `emmeans` (Estimated Marginal Means)

## Usage
1. Place the fitted circuit parameter data (`.csv`) into the `/data` folder.
2. Adjust the `custom_contrasts` vectors to match your specific factor levels.
3. Run `eis_statistical_analysis.R`. Outlier diagnostics and adjusted contrast tables will print directly to the console for review.

## Author
**Felipe Longaray Kadel**
*Materials Engineering Undergraduate - UFRGS*
* felipe.longaray@ufrgs.br





