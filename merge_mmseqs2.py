from Bio import SeqIO
import argparse
from pathlib import Path
from natsort import natsorted
import pickle
import os

def get_options():
    description = "Merges clusters from batched mmseqs2 runs."
    parser = argparse.ArgumentParser(description=description,
                                        prog='python merge_mmseqs.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--outpref',
                    default="clusters",
                    help='Output prefix. Default = "clusters"')
    IO.add_argument('--indir',
                    required=True,
                    help='Input directory for mmseqs2 results')
    IO.add_argument('--clusters',
                default=None,
                help='Cluster .pkl file generated from previous run of merge_mmseqs.py')
    return parser.parse_args()

def main():

    options = get_options()
    indir = options.indir
    outpref = options.outpref
    clusters = options.clusters

    # indir = "/hps/nobackup/jlees/mmseqs2/results"
    # outpref = "mmseqs2_batched_clusters"
    # clusters = "mmseqs2_batched_clusters.pkl"

    rep_files = []
    for path in Path(indir).glob("*_cluster.tsv"):
        # Print the path (file or directory) to the console
        rep_files.append(path)

    # sort based on batch number
    rep_files = natsorted(rep_files)

    if clusters is None:

        # fastas are hierarchically clustered, so need to extract the representatives in each file
        # iteratively and then determine which cluster they form in the next file
        rep_to_cluster = {}
        cluster_to_rep = {}
        for rep_file in rep_files:
            with open(rep_file, "r") as f:
                while True:
                    line = f.readline()
                    if not line:
                        break

                    split_line = line.rstrip().split("\t")
                    rep = split_line[0]
                    seq = split_line[1]

                    # if rep has not been seen before, add to new cluster
                    if rep not in rep_to_cluster:
                        rep_to_cluster[rep] = rep
                        cluster_to_rep[rep] = set()
                        cluster_to_rep[rep].add(rep)

                    # if sequence is in cluster_to_rep, means it has been 
                    # clustered with a new represenative
                    if seq in cluster_to_rep and seq != rep:

                        # update old reps
                        reps_set = cluster_to_rep[seq]

                        for prev_rep in reps_set:
                            rep_to_cluster[prev_rep] = rep

                        # account for duplicated IDs
                        try:
                            cluster_to_rep[rep].update(reps_set)
                            del cluster_to_rep[seq]
                        except KeyError:
                            pass
                    # elif seq != rep:
                    #     # add sequence to cluster
                    #     cluster_to_rep[rep].add(seq)
            print("Finished: {}".format(rep_file))

        print("No. clusters: {}".format(str(len(cluster_to_rep))))
        output_dicts = (rep_to_cluster, cluster_to_rep)

        with open(outpref + '.pkl', 'wb') as handle:
            pickle.dump(output_dicts, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # check results are correct
        print("Missing final centroids: ")
        with open(rep_files[-1], "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break

                split_line = line.rstrip().split("\t")
                rep = split_line[0]
                seq = split_line[1]

                if rep not in cluster_to_rep:
                    print(rep)
    
    else:
        # read in cluster file
        with open(clusters, 'rb') as handle:
            output_dicts = pickle.load(handle)
            rep_to_cluster, cluster_to_rep = output_dicts
        
        # write output
        with open(outpref + ".tsv", "w") as o:
            for rep_file in rep_files:
                current_rep = None
                with open(rep_file, "r") as f:
                    while True:
                        line = f.readline()
                        if not line:
                            break
                        split_line = line.rstrip().split("\t")

                        rep = split_line[0]
                        seq = split_line[1]

                        new_rep = rep_to_cluster[rep]

                        o.write(new_rep + "\t" + seq + "\n")

                print("Finished: {}".format(rep_file))


if __name__ == "__main__":
    main()