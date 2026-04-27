## training_preparation

### generate_embedding_labels.py

Generate classification labels for embedding visualization datasets. Creates standardized labeling schemes for training, validation, and novel genomes, supporting pretraining labels and PopPUNK cluster integration for downstream tasks.

**Usage:**
```
python generate_embedding_labels.py --original_labels <original.csv> --new_labels <new.csv> [--pretraining_labels <pretraining.csv>] [--clusters <clusters.csv>] [--outpref <prefix>]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--original_labels` | Yes | — | CSV file with genome names originally used in model training |
| `--new_labels` | Yes | — | CSV file with genome names not in model training |
| `--pretraining_labels` | No | None | CSV file with genome names used in model pretraining |
| `--clusters` | No | None | PopPUNK `clusters.csv` file, used if labels CSV has no second column |
| `--outpref` | No | `output` | Output prefix |

---

### match_files.py

Filter target files based on index file entries. Extracts matching filenames from target datasets using an index reference, supporting file path parsing and basename matching for data subset selection.

**Usage:**
```
python match_files.py --index <index_file> --target <target_file> [--outfile <output.txt>]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--index` | Yes | — | Path to the index file |
| `--target` | Yes | — | Path to the target file to search in |
| `--outfile` | No | `output.txt` | Output filename |

---

### merge_iter_poppunk.py

Merge results from multiple PopPUNK clustering iterations. Combines incremental clustering runs by mapping genomes to final cluster assignments, maintaining cluster consistency across iterative PopPUNK runs.

**Usage:**
```
python merge_iter_poppunk.py --infiles <run1.csv,run2.csv,...> [--outpref <prefix>]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--infiles` | Yes | — | Comma-separated list of PopPUNK run iteration files in order of running |
| `--outpref` | No | `merged_poppunk` | Output prefix |

---

### parse_sample_dates.py

Parse and stratify samples by collection dates. Processes ENA metadata files, supports date-based splitting by year, month, or day, handles incomplete date formatting, and enables temporal dataset organization.

**Usage:**
```
python parse_sample_dates.py --infile <ena_metadata> [--outpref <prefix>] [--downsample <genome_ids.txt>] [--split-by {year,month,day,none}] [--date-type {collection_date,last_updated}]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--infile` | Yes | — | ENA metadata file |
| `--outpref` | No | `output` | Output prefix |
| `--downsample` | No | None | File of genome IDs (no header) to downsample by |
| `--split-by` | No | `none` | How to stratify samples by date (`year`, `month`, `day`, or `none`) |
| `--date-type` | No | `collection_date` | Date field to use (`collection_date` or `last_updated`) |

---

### parse_xml.py

Parse NCBI BioSample XML files into tabular format. Extracts sample metadata including accession numbers, collection dates, geographic information, and taxonomic classifications from structured XML data.

**Usage:**
```
python parse_xml.py --infile <biosample.xml> [--outfile <output.tsv>]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--infile` | Yes | — | NCBI XML file from BioSample full XML download |
| `--outfile` | No | `output.tsv` | Output filename |

---

### remove_rare_tokens.py

Filter out rare tokens from genome datasets. Removes tokens occurring below specified frequency thresholds, supporting both absolute count and percentage-based filtering with options for incremental updates to existing filtered datasets.

**Usage:**
```
python remove_rare_tokens.py --infile <genomes_file> [--min <int>] [--min_perc <float>] [--counts <counts.pkl>] [--below_min_dict <dict.pkl>] [--outpref <prefix>]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--infile` | Yes | — | Input genome file |
| `--min` | No | `1` | Minimum token count to retain (tokens at or below this value are removed) |
| `--min_perc` | No | None | Minimum token frequency as a percentage (overrides `--min` if set) |
| `--counts` | No | None | Counts `.pkl` from a previous run (for incremental updates) |
| `--below_min_dict` | No | None | Below-min-dict `.pkl` from a previous run (for incremental updates) |
| `--outpref` | No | `genomes` | Output prefix |

---

### sample_PopPUNK.py

Sample representative genomes from PopPUNK clusters. Performs proportional or fixed-number sampling across clusters with minimum size constraints, ensuring balanced representation for training datasets.

**Usage:**
```
python sample_PopPUNK.py --clusters <clusters.csv> --genomes <poppunk_infile> [--prop <float>] [--number <float>] [--min <int>] [--outpref <prefix>] [--seed <int>]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--clusters` | Yes | — | Path to clusters file generated by PopPUNK |
| `--genomes` | Yes | — | Path to PopPUNK infile for generating sample paths |
| `--prop` | No | `0.1` | Proportion of each cluster to sample |
| `--number` | No | None | Fixed number of genomes to sample per cluster (overrides `--prop`) |
| `--min` | No | `100` | Minimum number of genomes to sample per cluster |
| `--outpref` | No | `sampled_genomes` | Output prefix |
| `--seed` | No | `42` | Random seed for reproducibility |

---

### stratify_N_genomes.py

Sample fixed number of genomes per cluster. Randomly selects specified number of genomes from each PopPUNK cluster while maintaining cluster assignments for stratified dataset creation.

**Usage:**
```
python stratify_N_genomes.py --clusters <clusters.csv> --genomes <tokenised_genomes> [--n-genomes <int>] [--outpref <prefix>] [--seed <int>]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--clusters` | Yes | — | Path to clusters file generated by PopPUNK |
| `--genomes` | Yes | — | Path to tokenised genomes file (genome per line with filename at start) |
| `--n-genomes` | No | `1` | Number of genomes to sample per cluster |
| `--outpref` | No | `stratified_genomes` | Output prefix |
| `--seed` | No | `42` | Random seed for reproducibility |

---

### stratify_data.py

Split genomic data into training, validation, and test sets. Performs lineage-aware stratified sampling to maintain population structure across data splits, with optional duplicate removal using distance matrices.

**Usage:**
```
python stratify_data.py --clusters <clusters.csv> --genomes <tokenised_genomes> [--val-size <float>] [--test-size <float>] [--outpref <prefix>] [--seed <int>] [--distances <dists_prefix>]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--clusters` | Yes | — | Path to clusters file generated by PopPUNK |
| `--genomes` | Yes | — | Path to tokenised genomes file (genome per line with filename at start) |
| `--val-size` | No | `0.1` | Proportion of genomes for validation |
| `--test-size` | No | `0.1` | Proportion of genomes for testing |
| `--outpref` | No | `stratified_genomes` | Output prefix |
| `--seed` | No | `42` | Random seed for reproducibility |
| `--distances` | No | None | Path to `.dists` prefix generated by PopPUNK (used to remove duplicates) |
