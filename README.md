# Latent Class Analysis of NHANES 2011–2018

Helen Hung · IDC 6940 · Advisor: Dr. Cohen

[Read the project](https://helenhung622.github.io/IDC6940_LCA/).

The report presents the primary analysis from LCA_MetS_v6: 5,507 adults,
seven binary indicators, one- through five-class model comparison, model
diagnostics, and descriptive characterization of the three-class solution.
The report does not include demographic predictor models or secondary
sample/indicator comparisons.

## Files

- `index.qmd`: report source, organized into Introduction, Methods, Analysis
  and Results, Conclusion, and Readings.
- `index.html`: rendered report served by GitHub Pages.
- `references.bib`: bibliography.
- `assets/`: two original v6 figures and three vector illustrations/charts.
- `tables/`: six aggregate v6 output tables, linked from the report.
- `analysis/lca_phase1.R`: Phase I subset of the v6 R script. It preserves
  the original model-fitting settings and clarifies figure legend labels.
- `REVISION_PLAN.md`: completed work and remaining checks (Traditional Chinese).

## Render the report

With Quarto installed, run from this repository:

```sh
quarto render index.qmd
```

The report uses the committed aggregate outputs and does not require
refitting the statistical models. Mathematical notation is rendered as
MathML, and images/styles are embedded in the generated HTML.

## Reproduce the Phase I analysis

The original cleaned participant-level CSV is not distributed in this
repository. Supply the same cleaned file and a separate output directory:

```sh
Rscript analysis/lca_phase1.R /absolute/path/to/nhanes2011_2018_adjw_clean_mets2.csv /absolute/path/to/output
```

Required R packages are declared at the start of the script. The script
writes tables, figures, and saved fit objects into the supplied output
directory. The published figures/tables are the existing v6 outputs; the
extracted script has not been rerun as part of this website revision.
Package versions and upstream preprocessing require further documentation.

Overall participant percentages use survey weights, strata and PSUs.
LCA and assigned-class summaries are unweighted. Class percentages are
modal-assignment percentages and should not be presented as national class
prevalence.

## Publishing

Render `index.qmd`, commit its HTML and any updated assets, tables,
bibliography, and style files, then push to the configured Pages branch.
The existing slide files are the original course template and have not
been updated in this revision.
