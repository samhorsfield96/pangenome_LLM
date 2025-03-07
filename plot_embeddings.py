
import argparse
import umap
import umap.plot
import pandas as pd
from collections import OrderedDict
import numpy as np
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def get_options():
    description = "Merges clusters from batched mmseqs2 runs."
    parser = argparse.ArgumentParser(description=description,
                                        prog='python merge_mmseqs.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--embeddings',
                    required=True,
                    help='Embeddings file generated by compute_sequence_embedding.py')
    IO.add_argument('--labels',
                    default=None,
                    help='csv file describing genome names in first column in same order as in embeddings file. No header. Can have second column with assigned clusters.')
    IO.add_argument('--clusters',
                    default=None,
                    help='PopPUNK clusters.csv file. Can be provided in addition to labels if no second column')
    IO.add_argument('--outpref',
                default="output",
                help='Output prefix. Default = "output"')
    return parser.parse_args()

def main():
    options = get_options()
    embeddings = options.embeddings
    labels = options.labels
    outpref = options.outpref
    clusters = options.clusters

    labels_dict = OrderedDict()
    if labels != None:
        print("Reading labels...")
        with open(labels, "r") as i1:
            for line in i1:
                split_line = line.rstrip().split(",")
                if len(split_line) >= 2:
                    labels_dict[split_line[0]] = split_line[1]
                else:
                    labels_dict[split_line[0]] = "NA"
    
        # only read clusters if labels also present
        if clusters != None:
            print("Reading clusters...")
            with open(clusters, "r") as i2:
                # read header
                i2.readline()
                for line in i2:
                    split_line = line.rstrip().split(",")
                    if split_line[0] in labels_dict:
                        labels_dict[split_line[0]] = split_line[1]

    
    # read embeddings
    print("Reading embeddings...")
    df = pd.read_csv(embeddings, header=None)
    #df.insert(loc=0, column='Cluster', value=cluster_list)
    #df.insert(loc=0, column='Sample', value=sample_list)
    
    reducer = umap.UMAP(random_state=42)

    print("Generating UMAP...")
    mapper = reducer.fit(df)
    UMAP_embedding = reducer.transform(df)

    cluster_list = [x for x in labels_dict.values()]
    sample_list = [x for x in labels_dict.keys()]
    #norm_cluster = np.array(cluster_list) / max(cluster_list)

    UMAP_embedding_df = pd.DataFrame(UMAP_embedding)
    UMAP_embedding_df.insert(loc=0, column='Cluster', value=cluster_list)
    UMAP_embedding_df.insert(loc=0, column='Sample', value=sample_list)

    UMAP_embedding_df.columns = ['Sample', 'Cluster', 'UMAP1', 'UMAP2']
    UMAP_embedding_df.to_csv(outpref + '.csv', index=False)

    print("Plotting UMAP...")
    if labels != None and cluster_list[0] != None:
        # unique_strings = list(set(cluster_list))
        # string_to_int = {s: i for i, s in enumerate(unique_strings)}
        # normalized_values = np.array([string_to_int[s] for s in cluster_list]) / (len(unique_strings) - 1)
        # cmap = cm.get_cmap("cubehelix")
        
        p = umap.plot.points(mapper, labels=UMAP_embedding_df['Cluster'], theme='fire')

        # plt.scatter(
        # UMAP_embedding[:, 0],
        # UMAP_embedding[:, 1],
        # c=normalized_values,
        # cmap=cmap)
        # plt.gca().set_aspect('equal', 'datalim')
        
    else:
        # plt.scatter(
        # UMAP_embedding[:, 0],
        # UMAP_embedding[:, 1])
        # plt.gca().set_aspect('equal', 'datalim')

        p = umap.plot.points(mapper)

    print("Saving file...")
    plt.savefig(outpref + ".png", dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    main()