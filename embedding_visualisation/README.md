## embedding_visualisation

### assign_strains_sparse.py

Assign strains using k-nearest neighbors classification on [ppsketchlib](https://github.com/bacpop/pp-sketchlib) distances. Performs strain classification by comparing query genomes to training set using sparse distance matrices, calculating accuracy metrics across multiple k values and generating detailed per-label performance statistics.

### generate_clusters_sparse.py

Generate Leiden clusters from sparse distance matrices using [ppsketchlib](https://github.com/bacpop/pp-sketchlib) distances. Performs community detection on genome similarity networks using various k-nearest neighbor and Leiden resolution parameters, optimizing clustering quality with Adjusted Rand Index and Adjusted Mutual Information scores.

### plot_embeddings.py

Plot UMAP visualizations of panBART embeddings. Creates dimensionality reduction plots from sequence embeddings with optional labeling by genome clusters, supporting color-coded visualization of population structure and configurable label limits for clarity.

### plot_embeddings_sparse.py

Plot UMAP visualizations from sparse [ppsketchlib](https://github.com/bacpop/pp-sketchlib) distance matrices. Generates 2D embeddings from sparse distance data with handling for isolated genomes, supports cluster labeling and creates publication-ready visualizations of genome relationships.
