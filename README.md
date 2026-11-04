# Salvage cryoablation after SBRT

Analysis code for a multicenter study of percutaneous cryoablation for in-field lung cancer recurrence after stereotactic body radiation therapy.

Patient data aren't included. To rerun the analysis, ask the study team for the data bundle and unzip it in this folder.

## Running it

You need R 4.5.0 and Python 3.12.

```sh
Rscript --vanilla reproduce.R
```

This installs the locked packages into the project folder, which needs internet access the first time, and writes results to `outputs/current/`.

To compare the results with the baseline in the bundle, install Poppler's `pdftoppm` and run:

```sh
Rscript --vanilla reproduce.R --verify
```

Numbers have to match on every platform. Figure files only have to match on macOS arm64, where the baseline was made and the only platform tested so far.

## Layout

| Path | Contents |
| --- | --- |
| `run_all.R` | The R analysis, in order |
| `analysis/R/` | Analysis modules |
| `analysis/scripts/` | Publication tables and the Python figures |
| `analysis/validation/` | Independent check in Python |
| `tests/` | Comparison with the baseline |

## Two sets of figures

The R figures in `publication_figures/` are the main set. The Python figures in `python_publication_figures/` were supplied separately and are reproduced as they were, so they differ slightly from the R analysis. Each run writes the differences to `statistical_consistency.csv`.
