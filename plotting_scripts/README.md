## plotting_scripts

### gene_freq_to_pseudolikelihood.R

Analyze relationship between gene frequency and pseudolikelihood scores. Merges gene frequency data with pseudolikelihood calculations, creates scatter plots with LOESS smoothing, and generates boxplots showing how gene frequency affects model predictions.

**Usage:**
```
Rscript gene_freq_to_pseudolikelihood.R
```

> **Note:** This script uses hardcoded input file paths. Edit the `filename` variables near the top of the script before running to point to the correct gene frequency and pseudolikelihood files.

---

### plot_loss.R

Visualize training metrics over epochs for model performance analysis. Creates multi-panel plots showing loss, precision, recall, and Cohen's kappa across training and validation datasets with publication-quality styling and logarithmic scaling for loss curves.

**Usage:**
```
Rscript plot_loss.R
```

> **Note:** This script uses a hardcoded input file path. Edit the `read.csv()` call near the top of the script to point to the parsed training log file produced by `parse_training_log.py`.

---

### plot_loss_model_selection.R

Compare model performance across different hyperparameter configurations. Generates bar plots of validation losses for parameter sets and uses random forest analysis to identify important hyperparameters for model selection.

**Usage:**
```
Rscript plot_loss_model_selection.R
```

> **Note:** This script uses a hardcoded input file path. Edit the `read.csv()` call near the top of the script to point to the grid search summary file (default: `synteny_grid_summary.txt`).

---

### plot_pseudolikelihood.R

Create comparative visualizations of pseudolikelihood scores across different conditions. Generates boxplots and statistical comparisons between trained, similar, and novel genome conditions with faceted plots for different inference types.

**Usage:**
```
Rscript plot_pseudolikelihood.R
```

> **Note:** This script uses a hardcoded input Excel file. Edit the `filename` variable near the top of the script to point to the pseudolikelihood results `.xlsx` file, and ensure the sheet names match your data.

---

### plot_token_likelihoods.R

Analyze token likelihood distributions across different lineages. Processes multiple pseudolikelihood files, calculates quantiles and summary statistics, and creates comparative boxplots with statistical significance testing between groups.

**Usage:**
```
Rscript plot_token_likelihoods.R
```

> **Note:** This script uses a hardcoded input directory. Edit the `indir` variable near the top of the script to point to the directory containing `*locus_pseudolikelihood.txt` files.
