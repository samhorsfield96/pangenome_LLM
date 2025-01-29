import argparse
from collections import defaultdict
import random
import os
import math
import numpy as np
import pickle

def parse_args():
    """
    Parse command-line arguments.

    This function parses the command-line arguments provided by the user and returns
    a Namespace object containing the parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Split data into training, validation and test datasets.")
    parser.add_argument("--clusters", type=str, required=True, help="Path to clusters file generated by PopPUNK.")
    parser.add_argument("--genomes", type=str, required=True, help="Path to tokenised genomes generated by tokenise_clusters.py. Must have genome per line with file name at beginning.")
    parser.add_argument("--val-size", type=float, default=0.1, help="Proportion of genomes for validation. Default = 0.1")
    parser.add_argument("--test-size", type=float, default=0.1, help="Proportion of genomes for testing. Default = 0.1.")
    parser.add_argument("--outpref", type=str, default="stratified_genomes", help="Output prefix. Default = stratified_genomes")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument('--distances', help='Path to .dists prefix generated by PopPUNK (required)',
                                required=True)

    args = parser.parse_args()

    return args

def main():

    args = parse_args()
    clusters = args.clusters
    genomes = args.genomes
    val_size = args.val_size
    test_size = args.test_size
    outpref = args.outpref
    seed = args.seed
    train_size = 1.0 - (val_size + test_size)
    distances = args.distances

    # set seed
    random.seed(seed)

    if train_size <= 0.0:
        print("Sum of --val-size and --test-size cannot be more than 1.0")
        sys.exit(1)

    cluster_dict = defaultdict(list)
    with open(clusters, "r") as f:
        f.readline()
        for line in f:
            split_line = line.strip().split(",")
            
            # add genome to cluster entry 
            cluster_dict[split_line[1]].append(split_line[0])

    # open stored distances
    with open(distances + ".pkl", 'rb') as pickle_file:
        rlist, qlist, self = pickle.load(pickle_file)
        r_names = [os.path.splitext(os.path.basename(name))[0] for name in rlist]
        #q_names = [os.path.splitext(os.path.basename(name))[0] for name in qlist]

    # Load matrix
    X = np.load(distances + ".npy")
    zero_dists = np.where(X.sum(axis=1) == 0.0)[0]
    del X
    #print(zero_dists)

    # identify indices of zero distances
    zero_dist_indices = []
    num_genomes = len(r_names)
    duplicate_set = set()
    #duplicate_dict = defaultdict(set)
    rows, cols = np.triu_indices(num_genomes, k=1)
    for dist in zero_dists:
        i, j = rows[dist], cols[dist]
        duplicate_set.add(r_names[j])
        #duplicate_dict[r_names[i]].add(r_names[j])

    #print(duplicate_dict) 

    # for each cluster, shuffle and randomly assign to training, don't keep tally as just all remaining genomes are training
    cluster_len_dict = {}
    for cluster in cluster_dict.keys():
        genome_list = cluster_dict[cluster]
        #print_on = False

        # go through genome list and remove duplicates
        #dup_set = set()
        #for genome in genome_list:
            #if genome in duplicate_set:
                #print_on = True
                #dup_set.update(duplicate_dict[genome])
                #print(f"{genome} : {duplicate_dict[genome]}")
        
        #if print_on:
            #print("pre")
            #print(dup_set)
            #print(genome_list)
        genome_list = [genome for genome in genome_list if genome not in duplicate_set]

        #if print_on:
            #print("post")
            #print(genome_list)
        
        # keep track of original cluster size 
        cluster_len_dict[cluster] = len(genome_list)

        # calculate number of genomes to move and then shuffle
        num_to_move = math.ceil(len(genome_list) * train_size)
        random.shuffle(genome_list)

        # remove genomes to new list
        cluster_dict[cluster] = genome_list[num_to_move:]

    # for each cluster, shuffle and randomly assign to validation
    val_genomes = set()
    for cluster in cluster_dict.keys():
        genome_list = cluster_dict[cluster]

        # calculate number of genomes to move and then shuffle
        num_to_move = math.ceil(cluster_len_dict[cluster] * val_size)
        random.shuffle(genome_list)

        # move genomes to new list
        val_genomes.update(genome_list[:num_to_move])
        cluster_dict[cluster] = genome_list[num_to_move:]
    
    # for remaining genomes are for testing
    test_genomes = set()
    for cluster, genome_list in cluster_dict.items():
        test_genomes.update(genome_list)

    # sample genomes from input file
    with open(outpref + "_train_genomes.txt", 'w') as o_train, open(outpref + "_val_genomes.txt", 'w') as o_val, open(outpref + "_test_genomes.txt", 'w') as o_test:
        with open(genomes, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                
                genome_id = os.path.splitext(line.rstrip().split("\t")[0])[0]

                if genome_id in test_genomes:
                    o_test.write(line)
                elif genome_id in val_genomes:
                    o_val.write(line)
                else:
                    o_train.write(line)
                
if __name__ == "__main__":
    main()