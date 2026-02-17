import argparse
from collections import defaultdict

def get_options():
    description = "Merge results from multiple PopPUNK run iterations."
    parser = argparse.ArgumentParser(description=description,
                                        prog='python merge_iter_poppunk.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--infiles',
                    required=True,
                    help='Comma separated list of PopPUNK run iterations in order of running.')
    IO.add_argument('--outpref',
                default="merged_poppunk",
                help='Output prefix.')

    return parser.parse_args()

def get_cluster_set(infile):
    clusters_set = set()
    with open(infile, "r") as f:
        header = f.readline()
        while True:
            line = f.readline()
            if not line:
                break

            cluster = line.rstrip().split(",")[1]
            clusters_set.add(cluster)
    return clusters_set

def main():
    options = get_options()
    infiles = options.infiles.split(",")
    outpref = options.outpref
    
    # get clusters in last file
    final_clusters_dict = {x: set(x.split("_")) for x in get_cluster_set(infiles[-1])}

    # iterate backwards over each infile to assign genomes to the same cluster, as final file is most clustered
    cluster_dict = defaultdict(set)
    file_count = 1
    for infile in infiles:
        # generate a set from the cluster IDs for each cluster unique cluster, then iterate through each line, grouping genomes together. 
        # Then determine whether the clusters have any itersection with those cluster names. 
        # If they do, assignment them to the cluster, and work backwards up the iteration.
        with open(infile, "r") as f:
            header = f.readline()
            while True:
                line = f.readline()
                if not line:
                    break

                split_line = line.rstrip().split(",")
                genome = split_line[0]
                cluster = split_line[1]

                cluster_set = set(cluster.split("_"))

                for key, entry in final_clusters_dict.items():
                    if len(entry.intersection(cluster_set)) > 0:
                        cluster_dict[key].add(genome)
                        break
        print(f"Finished file {file_count}")
        file_count += 1

    # write output
    with open(outpref + ".csv", "w") as f:
        f.write("Taxon,Cluster\n")
        for key, entry in cluster_dict.items():
            for genome in entry:
                f.write(f"{genome},{key}\n")


if __name__ == "__main__":
    main()