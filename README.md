# Loblaw Bio Immune Cell Analysis

This project loads immune cell count data into SQLite, computes relative cell
population frequencies, compares melanoma miraclib responders with
non-responders, and summarizes a baseline melanoma PBMC subset.

## Dashboard

After starting the dashboard, open:

[http://localhost:8000](http://localhost:8000)

In GitHub Codespaces, use the forwarded port URL for port `8000`.

## Reproduce the Pipeline in GitHub Codespaces

The project uses only Python standard library modules.

```bash
make setup
make pipeline
make dashboard
```

`make pipeline` runs the full workflow from start to finish:

1. Creates and loads `cell_counts.db` from `cell-count.csv`.
2. Generates per-sample cell population frequencies.
3. Runs responder vs non-responder statistics for melanoma PBMC samples from
   miraclib-treated patients.
4. Generates the baseline melanoma PBMC miraclib subset summaries.
5. Builds the HTML dashboard in `dashboard/index.html`.

## Input and Output Files

Input:

- `cell-count.csv`: source data with sample metadata and counts for `b_cell`,
  `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, and `monocyte`.

Generated outputs:

- `cell_counts.db`: SQLite database.
- `cell_population_frequencies.csv`: Part 2 relative frequency table.
- `miraclib_melanoma_pbmc_stats.csv`: Part 3 statistical results.
- `miraclib_melanoma_pbmc_boxplot.svg`: Part 3 responder/non-responder boxplot.
- `baseline_melanoma_miraclib_pbmc_samples.csv`: Part 4 baseline sample list.
- `baseline_project_summary.csv`: Part 4 samples per project.
- `baseline_response_summary.csv`: Part 4 responder/non-responder subject counts.
- `baseline_gender_summary.csv`: Part 4 male/female subject counts.
- `dashboard/index.html`: local dashboard.

## Relational Database Schema

The SQLite schema is created by `load_data.py` and uses four tables:

- `patients`: one row per subject. Stores project, indication, age, gender,
  treatment, and response.
- `samples`: one row per biological sample. Stores sample type and
  `time_from_treatment_start`, with a foreign key to `patients`.
- `cell_populations`: one row per immune population name.
- `cell_counts`: one row per sample/population pair, storing the observed cell
  count.

This design separates subject metadata, sample metadata, population definitions,
and measurements. It avoids repeating all subject metadata for each cell
population count, and it makes new analytics easy to express with joins. For
example, Part 3 filters on patient and sample attributes, then joins to
`cell_counts` for measurements.

The schema also scales well if the dataset grows to hundreds of projects and
thousands or millions of samples. Adding more projects only adds rows to
`patients`; adding more samples adds rows to `samples`; adding new immune cell
types adds rows to `cell_populations`; and measurements continue to append to
`cell_counts`. Indexes on sample subject IDs, treatment time, and population IDs
support common filters and joins. If analytics expanded further, this model
could support additional indexes, database views for reusable cohorts, or
separate tables for derived features without changing the raw measurement
structure.

## Code Structure

- `load_data.py`: initializes the SQLite schema and loads all rows from
  `cell-count.csv`.
- `analyze_data.py`: computes relative frequency percentages for each
  sample/population pair and writes `cell_population_frequencies.csv`.
- `statistical_analysis.py`: filters melanoma PBMC samples from miraclib-treated
  patients, compares responders and non-responders using Mann-Whitney U tests,
  applies Benjamini-Hochberg FDR correction, and writes the statistics CSV and
  SVG boxplot.
- `subset_analysis.py`: identifies baseline melanoma PBMC samples from
  miraclib-treated patients and writes the requested project, response, and
  gender summaries.
- `build_dashboard.py`: builds a static HTML dashboard from the generated CSV,
  SVG, and database outputs.
- `Makefile`: provides the required `setup`, `pipeline`, and `dashboard`
  targets for automated grading.

The scripts are intentionally small and sequential. Each script owns one
analysis step, writes durable outputs, and reads from the SQLite database or
previous pipeline outputs. This makes the pipeline easy to grade, rerun, and
debug in Codespaces.

## Statistical Summary

For melanoma PBMC samples from miraclib-treated patients, no immune cell
population was significant after Benjamini-Hochberg FDR correction at 0.05.
`cd4_t_cell` had the strongest unadjusted signal, with higher median relative
frequency in responders, but its FDR-adjusted p-value was above 0.05.
