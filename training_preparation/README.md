## training_preparation

### generate_embedding_labels.py

Generate classification labels for embedding visualization datasets. Creates standardized labeling schemes for training, validation, and novel genomes, supporting pretraining labels and PopPUNK cluster integration for downstream tasks.

### match_files.py

Filter target files based on index file entries. Extracts matching filenames from target datasets using an index reference, supporting file path parsing and basename matching for data subset selection.

### merge_iter_poppunk.py

Merge results from multiple PopPUNK clustering iterations. Combines incremental clustering runs by mapping genomes to final cluster assignments, maintaining cluster consistency across iterative PopPUNK runs.

### parse_sample_dates.py

Parse and stratify samples by collection dates. Processes ENA metadata files, supports date-based splitting by year, month, or day, handles incomplete date formatting, and enables temporal dataset organization.

### parse_xml.py

Parse NCBI BioSample XML files into tabular format. Extracts sample metadata including accession numbers, collection dates, geographic information, and taxonomic classifications from structured XML data.

### remove_rare_tokens.py

Filter out rare tokens from genome datasets. Removes tokens occurring below specified frequency thresholds, supporting both absolute count and percentage-based filtering with options for incremental updates to existing filtered datasets.

### sample_PopPUNK.py

Sample representative genomes from PopPUNK clusters. Performs proportional or fixed-number sampling across clusters with minimum size constraints, ensuring balanced representation for training datasets.

### stratify_N_genomes.py

Sample fixed number of genomes per cluster. Randomly selects specified number of genomes from each PopPUNK cluster while maintaining cluster assignments for stratified dataset creation.

### stratify_data.py

Split genomic data into training, validation, and test sets. Performs lineage-aware stratified sampling to maintain population structure across data splits, with optional duplicate removal using distance matrices.
