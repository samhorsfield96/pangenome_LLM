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
                    default="clusters.tsv",
                    help='Output filename. Default = "clusters.tsv"')
    IO.add_argument('--indir',
                    required=True,
                    help='Input directory for mmseqs2 results')
    IO.add_argument('--final-clusters',
                required=True,
                help='Merged representatives from running MMseqs2 on all representatives from all files in indir.')
    IO.add_argument('--clusters',
                default=None,
                help='Cluster .pkl file generated from previous run of merge_mmseqs.py')
    return parser.parse_args()

def merge_mmseqs2(rep_files, final_cluster, clusters):
    # sort based on batch number
    rep_files = natsorted(rep_files)

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
                    print(f"Error: {seq} already in cluster_to_rep.")
                    print(f"cluster_to_rep[{seq}]: ")
                    print(cluster_to_rep[seq])

        print("Finished: {}".format(rep_file))

    # Now go through final representatives and reassign clusters
    with open(final_cluster, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break

            split_line = line.rstrip().split("\t")
            rep = split_line[0]
            seq = split_line[1]

            # if rep has not been seen before, add to new cluster
            if rep not in rep_to_cluster:
                print(f"Error: final {rep} not in previous cluster.")

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

    print("No. clusters: {}".format(str(len(cluster_to_rep))))
    output_dicts = (rep_to_cluster, cluster_to_rep)

    with open(clusters, 'wb') as handle:
        pickle.dump(output_dicts, handle, protocol=pickle.HIGHEST_PROTOCOL)

def write_mmseqs2_clusters(rep_files, outfile, clusters):
    # sort based on batch number
    rep_files = natsorted(rep_files)

    # read in cluster file
    with open(clusters, 'rb') as handle:
        output_dicts = pickle.load(handle)
        rep_to_cluster, cluster_to_rep = output_dicts
    
    # write output
    with open(outfile, "w") as o:
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

def main():

    options = get_options()
    indir = options.indir
    outfile = options.outfile
    final_clusters = options.final_clusters
    clusters = options.clusters

    # indir = "/hps/nobackup/jlees/mmseqs2/results"
    # outpref = "mmseqs2_batched_clusters"
    # clusters = "mmseqs2_batched_clusters.pkl"

    rep_files = []
    for path in Path(indir).glob("*_cluster.tsv"):
        # Print the path (file or directory) to the console
        rep_files.append(path)

    # sort based on batch number, remove final file
    rep_files = [file for file in rep_files if file != final_clusters]
    rep_files = natsorted(rep_files)

    if clusters == None:
        merge_mmseqs(rep_files, final_clusters, outfile)
    else:
        write_mmseqs2_clusters(rep_files, outfile, clusters)

if __name__ == "__main__":
    main()