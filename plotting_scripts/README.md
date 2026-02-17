## plotting_scripts

### gene_freq_to_pseudolikelihood.R

Analyze relationship between gene frequency and pseudolikelihood scores. Merges gene frequency data with pseudolikelihood calculations, creates scatter plots with LOESS smoothing, and generates boxplots showing how gene frequency affects model predictions.

### plot_loss.R

Visualize training metrics over epochs for model performance analysis. Creates multi-panel plots showing loss, precision, recall, and Cohen's kappa across training and validation datasets with publication-quality styling and logarithmic scaling for loss curves.

### plot_loss_model_selection.R

Compare model performance across different hyperparameter configurations. Generates bar plots of validation losses for parameter sets and uses random forest analysis to identify important hyperparameters for model selection.

### plot_pseudolikelihood.R

Create comparative visualizations of pseudolikelihood scores across different conditions. Generates boxplots and statistical comparisons between trained, similar, and novel genome conditions with faceted plots for different inference types.

### plot_token_likelihoods.R

Analyze token likelihood distributions across different lineages. Processes multiple pseudolikelihood files, calculates quantiles and summary statistics, and creates comparative boxplots with statistical significance testing between groups.
