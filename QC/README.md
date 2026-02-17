## QC

### LLM_accuracy.py

Compares generated genome sequences to prompted reference genomes. Calculates identity between generated and reference sequences using pairwise alignment, handling reverse complement sequences and providing detailed accuracy metrics per gene token.

### calc_gene_frequency.py

Calculate token frequency across the population. Computes gene presence/absence frequencies, handles paralog counting, and can combine forward/reverse strand genes using absolute IDs. Optionally integrates with annotation data for functional frequency analysis.

### compare_mmseqs2_clusters.py

Compare clustering results between different MMseqs2 runs. Analyzes identical, partial, and non-matching clusters between query and database cluster files, with options to filter by minimum cluster size and cache intermediate results for efficiency.

### parse_training_log.py

Parse training logs from machine learning model training. Extracts training, validation, and test metrics including loss, perplexity, accuracy, precision, recall, and F1 scores, organizing them into structured format for analysis and visualization.
