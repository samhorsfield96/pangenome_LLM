## sequence_generation

### generate_genomes.py

Generate synthetic genomes from PanBART models. Integrates synteny predictions with gene-level sequence generation to create complete genome assemblies with configurable genome sizes and sampling parameters.

**Usage:**
```
python generate_genomes.py --reps <reps.pkl> --syntenyLLM <synteny_model> --synteny_tokeniser <synteny.bin> --geneLLM <gene_model> --gene_tokeniser <genes.bin> --nanoGPT <nanoGPT_dir> [--outdir <dir>] [--num_samples <int>] [--max_genome_size <int>] [--min_genome_size <int>] [--max_gene_prop <float>] [--temperature <float>] [--top_k <int>] [--seed <int>] [--device <device>] [--no_compile]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--reps` | Yes | — | Output `.pkl` from `group_sequences.py` |
| `--syntenyLLM` | Yes | — | Path to nanoGPT model trained on gene synteny |
| `--synteny_tokeniser` | Yes | — | Path to synteny tokeniser `.bin` |
| `--geneLLM` | Yes | — | Path to nanoGPT model trained on gene sequences |
| `--gene_tokeniser` | Yes | — | Path to gene tokens `.bin` |
| `--nanoGPT` | Yes | — | Path to nanoGPT directory |
| `--outdir` | No | `predictions` | Output directory |
| `--num_samples` | No | `1` | Number of genomes to generate |
| `--max_genome_size` | No | `3000` | Maximum number of gene tokens to generate |
| `--min_genome_size` | No | `100` | Minimum number of gene tokens to generate |
| `--max_gene_prop` | No | `1.5` | Maximum length proportion of simulated gene relative to representative |
| `--temperature` | No | `0.8` | Sampling randomness (1.0 = unchanged, <1.0 = less random, >1.0 = more random) |
| `--top_k` | No | `200` | Retain only the top-k most likely tokens |
| `--seed` | No | None | Random seed for sequence sampling |
| `--device` | No | `cuda` | Device to use (e.g. `cpu`, `cuda`, `cuda:0`) |
| `--no_compile` | No | False | Disable PyTorch 2 model compilation |

---

### synteny_accuracy.py

Evaluate synteny accuracy of generated genomes using network analysis. Compares gene neighborhood relationships between generated and reference genomes using graph-based metrics to assess preservation of genomic context and gene order.

**Usage:**
```
python synteny_accuracy.py --real <real_genomes.txt> --sim <simulated_genomes.txt> [--outpref <prefix>]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--real` | Yes | — | Path to real genome data `.txt` file |
| `--sim` | Yes | — | Path to simulated genome data `.txt` file |
| `--outpref` | No | `comparisons` | Output prefix |
