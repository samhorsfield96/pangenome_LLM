import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Get positions of specific tokens.")
    parser.add_argument('--infile', required=True, help='Path to tokenised genome data .txt file.')
    parser.add_argument("--tokens", type=str, required=True, help="Comma separated list of tokens to search for.")
    parser.add_argument("--outpref", default="plot", help="Output prefix")
    parser.add_argument('--absolute-IDs',
                    action="store_true",
                    default=False,
                    help='Use absolute gene IDs to calculate frequency, combining forward and reverse strand genes.')

    args = parser.parse_args()

    return args

def main():

    args = parse_args()
    infile = args.infile
    tokens = args.tokens.split(",")
    outpref = args.outpref
    absolute_IDs = args.absolute_IDs

    genome_token_dict = {}
    with open(infile, "r") as f1:
        while True:
            line = f1.readline()
            if not line:
                break
            gene_presence_genome = {}

            split_line = line.rstrip().split("\t")
            genome_id = split_line[0]

            split_line2 = split_line[-1].split(" ")
            token_dict = {token: [] for token in tokens}
            for idx, gene in enumerate(split_line2):
                if gene != "_":
                    abs_gene_ID = str(abs(int(gene)))
                    gene = abs_gene_ID if absolute_IDs else gene
                
                if gene in token_dict:
                    token_dict[gene].append(idx)

            genome_token_dict[genome_id] = token_dict


    with open(outpref + "_token_locations.csv", "w") as f:
        for genome_id, token_dict in genome_token_dict.items():
            for gene, idx_list in token_dict.items():
                for idx in idx_list:
                    f.write(f"{genome_id},{gene},{idx}\n")
                
    
if __name__ == "__main__":
    main()