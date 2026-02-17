## token_annotation

### annotate_tokens.py

Annotate tokenized gene clusters with functional information. Maps gene cluster representatives to annotation data, creating comprehensive token annotation files that link numerical tokens to gene names and functional descriptions.

### annotation_to_token.py

Match tokenizer vocabulary with annotated token data. Filters annotation files to include only tokens present in the tokenizer vocabulary, ensuring compatibility between model inputs and functional annotations.

### find_concurrent_genes.py

Identify genomes containing specific gene combinations. Searches for co-occurring genes within genomes, supporting single gene searches or pairwise gene combinations, and outputs filtered genome subsets for downstream analysis.

### get_gene_positions.py

Extract positional information for specific tokens across genomes. Maps token locations within genome sequences, supporting absolute ID conversion and generating position matrices for spatial analysis of gene distributions.

### token_likelihood_position.py

Analyze token likelihood distributions across genomic positions. Calculates position-specific token frequencies, generates statistical profiles, and creates visualizations of token positional patterns within genome assemblies.
