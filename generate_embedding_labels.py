
import argparse
from collections import OrderedDict

def get_options():
    description = "Generates labels for embedding plots."
    parser = argparse.ArgumentParser(description=description,
                                        prog='python generate_embedding_labels.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--original_labels',
                    required=True,
                    help='csv file describing genome names orignally in model training in first column.')
    IO.add_argument('--new_labels',
                    required=True,
                    help='csv file describing genome names not in model training in first column.')
    IO.add_argument('--clusters',
                    default=None,
                    help='PopPUNK clusters.csv file. Can be provided in addition to labels if no second column')
    IO.add_argument('--outpref',
                default="output",
                help='Output prefix. Default = "output"')
    return parser.parse_args()

def main():
    options = get_options()
    original_labels = options.original_labels
    new_labels = options.new_labels
    outpref = options.outpref
    clusters = options.clusters

    # key for genome identification
    # 1 = training genome
    # 2 = new genome matching strain found in training
    # 3 = new genome not matching training strain

    labels_dict = OrderedDict()
    original_clusters = set()
    known_isolates = set()
    print("Reading original labels...")
    with open(original_labels, "r") as i1:
        for line in i1:
            split_line = line.rstrip().split(",")
            if len(split_line) >= 2:
                labels_dict[split_line[0]] = (split_line[1], 1)
                original_clusters.add(split_line[1])
            else:
                labels_dict[split_line[0]] = (None, 1)
    
        # only read clusters if labels also present
        if clusters != None:
            print("Reading clusters...")
            with open(clusters, "r") as i2:
                # read header
                i2.readline()
                for line in i2:
                    split_line = line.rstrip().split(",")
                    if split_line[0] in labels_dict:
                        labels_dict[split_line[0]] = (split_line[1], 1)
                        original_clusters.add(split_line[1])
                        known_isolates.add(split_line[0])

    print("Reading new labels...")
    new_genomes = set()
    with open(new_labels, "r") as i1:
        for line in i1:
            split_line = line.rstrip().split(",")

            # if already observed, ignore
            if split_line[0] not in labels_dict:
                if len(split_line) >= 2:
                    if split_line[1] in original_clusters:
                        labels_dict[split_line[0]] = (split_line[1], 2)
                    else:
                        labels_dict[split_line[0]] = (split_line[1], 3)
                else:
                    labels_dict[split_line[0]] = (None, None)
                new_genomes.add(split_line[0])
    
        # only read clusters if labels also present
        if clusters != None:
            print("Reading clusters...")
            with open(clusters, "r") as i2:
                # read header
                i2.readline()
                for line in i2:
                    split_line = line.rstrip().split(",")
                    # only update if new genome
                    if split_line[0] in new_genomes:
                        known_isolates.add(split_line[0])
                        if split_line[1] in original_clusters:
                            labels_dict[split_line[0]] = (split_line[1], 2)
                        else:
                            labels_dict[split_line[0]] = (split_line[1], 3)
    
    # remove isolates not part of dataset
    if len(known_isolates) > 0:
        to_remove = set()
        for isolate_ID, cluster_set in labels_dict.items():
            if isolate_ID not in known_isolates:
                to_remove.add(isolate_ID) 
        
        for isolate_ID in to_remove:
            del labels_dict[isolate_ID]

    with open(outpref + ".csv", "w") as o:
        for isolate_ID, cluster_set in labels_dict.items():
            o.write(str(isolate_ID) + "," + str(cluster_set[1]) + "," + str(cluster_set[0]) + "\n")

if __name__ == "__main__":
    main()