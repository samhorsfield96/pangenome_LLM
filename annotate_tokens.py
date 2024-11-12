from Bio import SeqIO
import argparse
import pickle

def get_options():
    description = "Annotates tokenised gene clusters"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python annotate_tokens.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--clusters',
                    required=True,
                    help='Cluster file.')
    IO.add_argument('--annotations',
            required=True,
            help='Path to gene annotation FASTA file.')
    IO.add_argument('--outpref',
                default="annotated_tokens",
                help='Output prefix.')
    return parser.parse_args()

def main():
    options = get_options()
    clusters = options.clusters
    annotations = options.annotations
    outpref = options.outpref
    
    annotation_dict = {}
    # read in annotations
    fasta_sequences = SeqIO.parse(open(annotations),'fasta')
    for fasta in fasta_sequences:
        name, description = fasta.id, fasta.description
        
        # parse name
        parsed_name = name[3:].replace(".contig", "_")
        
        annotation_dict[parsed_name] = description
        #print(parsed_name)

    # read in cluster file
    with open(clusters, 'rb') as handle:
        reps_dict = pickle.load(handle)
    
    # iterate through reps_dict, assigning details to each entry
    with open(outpref + ".tsv", "w") as o:
        o.write("Name\tToken\tAnnotation\n")
        for token, name in reps_dict.items():
            annotation = annotation_dict[name]

            o.write("{}\t{}\t{}\n".format(str(name), str(token), str(annotation)))


if __name__ == "__main__":
    main()