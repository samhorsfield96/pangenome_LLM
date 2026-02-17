from scipy.spatial.distance import pdist
from scipy.spatial import distance
import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import networkx as nx
from scipy.cluster.hierarchy import single, fcluster, dendrogram



def main():
    #parse options
    #options = get_options()
    #infile = options.infile
    #plot = options.plot
    #outpref = options.outpref
    #cutoff = options.cutoff

    #for debugging
    infile = "tokenised_genomes.txt"
    plot = True
    outpref = "synteny"
    subsample = 100000
    kde_bandwidth = 0.01
    kde_tolerance = 0.2

    genome_list = []
    genome_name_list = []
    max_value = 0
    n_genomes = 0
    print("Reading inputs...")
    with open(infile, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            split_genome = line.rstrip().split("\t")
            genome_name_list.append(split_genome[0])
            genome = split_genome[1].split(" ")
            
            # ensure gene ids are absolute numbers, make 0 indexed
            genome = set([abs(int(x)) - 1 for x in genome if x != "_"])
            genome_list.append(genome)

            # determine maximum column number
            max_gene = max(genome) + 1
            if max_gene > max_value:
                max_value = max_gene
            
            n_genomes += 1

    # generate numpy array where each row is genome and each column is gene
    pangenome_array = np.zeros((n_genomes, max_value), dtype=int)

    for genome_idx, genome in enumerate(genome_list):
        genome_array = np.zeros((max_value), dtype=int)
        genome_array[list(genome)] = 1
        pangenome_array[genome_idx] = genome_array

    print("Calculating Jaccard distances...")
    jaccard_acc = pdist(pangenome_array, metric='jaccard')
    #print(len(jaccard_acc))

    if plot == True:
        # Plot the data and density function with identified peaks
        plt.figure(figsize=(12, 6))
        plt.hist(jaccard_acc, bins=50, density=True, alpha=0.5, label='Histogram of Data')

        plt.title('Data Distribution')
        plt.xlabel('Accessory Distance')
        plt.ylabel('Density')
        plt.legend()
        plt.show()

        plt.savefig(outpref + "_hist.png")

    # optimise bandiwth for Kernel density
    # subsample to N
    jaccard_acc_sample = np.random.choice(jaccard_acc, subsample)

    
    print("Fitting KDE model...")
    # bandwidth grid search
    #bandwidths = 10 ** np.linspace(-2, 1, 10)
    #grid = GridSearchCV(KernelDensity(kernel='gaussian'),{'bandwidth': bandwidths})
    #grid.fit(jaccard_acc_sample[:, None])
    #kde = grid.best_estimator_
    
    #print(kde)
    kde = KernelDensity(bandwidth=kde_bandwidth, kernel='epanechnikov')
    kde.fit(jaccard_acc_sample[:, None])
    x = np.linspace(0, 1000, 1000)[:, None] / 1000
    log_density = kde.score_samples(x)
    density = np.exp(log_density)

    # Find local maxima
    maxima_indices = argrelextrema(density, np.greater)
    cluster_centers = x[maxima_indices]
    cluster_centers = cluster_centers.T[0]
    cluster_centers = np.sort(cluster_centers)

    # get first and second cluster
    print("Assigning clusters...")
    clusters = np.zeros(len(jaccard_acc), dtype=int)
    print("Cluster centers: {}".format(cluster_centers))
    if len(cluster_centers) > 1:
        # decide which clusters to compare
        cluster1 = cluster_centers[0]

        # ensure that next peak is sufficiently different from first peak
        for i in range(1, len(cluster_centers)):
            cluster2 = cluster_centers[i]
            if cluster2 >= cluster1 * (1 + kde_tolerance):
                break

        cutoff = (cluster1 + cluster2) / 2
        print("Accessory Jaccard cutoff: {}".format(cutoff))
        clusters = np.zeros(len(jaccard_acc), dtype=int)
        # assign as cluster 1 if above cutoff
        mask = jaccard_acc > cutoff
        clusters[mask] = 1
    else:
        print("No clusters detected. All same strain.")

    #print(cluster_centers)

    # assign clusters
    #clusters = np.abs(jaccard_acc[:, np.newaxis] - cluster_centers).argmin(axis=1)

    # generate square plot to allow indexing
    clusters = distance.squareform(clusters)
    #print(clusters)
    
    if plot == True:
        # Plot the data and density function with identified peaks
        plt.close()
        plt.figure(figsize=(12, 6))
        plt.hist(jaccard_acc, bins=50, density=True, alpha=0.5, label='Histogram of Data')
        plt.plot(x, density, label='Kernel Density Estimation')
        plt.scatter(x[maxima_indices], density[maxima_indices], color='red', zorder=5, label='Peaks')

        plt.title('Kernel Density Estimation and Data Distribution')
        plt.xlabel('Accessory Distance')
        plt.ylabel('Density')
        plt.legend()
        plt.show()

        plt.savefig(outpref + "_kernel_plot.png")

    #print(clusters)

    # iterate over upper half with diagonal
    rows, columns = clusters.shape
    #print(rows)
    #print(columns)
    G = nx.Graph()

    print("Generating network...")

    for row in range(rows):
        for col in range(row, columns):
            # only need to detect 0s as this is same cluster
            if clusters[row, col] == 0:
                genome1 = genome_name_list[row]
                genome2 = genome_name_list[col]

                G.add_edge(genome1, genome2)

    graph_clusters = list(nx.connected_components(G))

    print("Writing output...")
    #print(graph_clusters)
    # write out clusters of genomes
    with open(outpref + "_clusters.txt", "w") as f:
        for idx, cluster in enumerate(graph_clusters):
            for genome in cluster:
                f.write(str(idx) + "\t" + genome + "\n")

if __name__ == "__main__":
    main()