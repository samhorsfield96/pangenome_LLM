import argparse
import json

def get_options():
    description = "Get the disribution of counts of a given token in the genomes."
    parser = argparse.ArgumentParser(description=description,
                                        prog='python token_count_distribution.py')
    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--genomes',
                    required=True,
                    help='Tokenised input genomes from WTBcluster')
    IO.add_argument('--tokeniser',
                    required=True,
                    help='Tokeniser from panBART fitting.')
    IO.add_argument('--tokens',
                    default="",
                    help='Comma separated list of tokens to generate distribution for. If not specified, just looks for unknowns and contig breaks.')
    IO.add_argument('--outpref',
                default="output",
                help='Output prefix. Default = "output"')
    return parser.parse_args()


def main():
    options = get_options()

    genomes = options.genomes
    tokeniser = options.tokeniser
    tokens = set((options.tokens + ",<unk>,_").split(","))
    tokens = [str(x) for x in tokens]
    outpref = options.outpref

    with open(tokeniser) as f:
        token_dict = json.load(f)["model"]["vocab"]

    token_str_to_token = {}
    for token in tokens:
        if str(token) in token_dict:
            token_str_to_token[str(token)] = token_dict[str(token)]

    # read through files to identify counts of tokens
    token_count_list = {x: [] for x in token_str_to_token.keys()}
    with open(genomes, "r") as f1:
        while True:
            line = f1.readline()
            if not line:
                break

            # hold counts in dictionary
            token_counts = {x: 0 for x in token_str_to_token.keys()}

            stripped_line = line.rstrip().split("\t")
            genome_id = stripped_line[0]
            genome_entry = stripped_line[-1]

            split_line = genome_entry.split(" ")
            for gene in split_line:
                if gene in token_dict and gene in token_counts:
                    token_counts[gene] += 1
                elif gene not in token_dict:
                    token_counts["<unk>"] += 1

            for key, entry in token_counts.items():
                token_count_list[key].append(entry)

    with open(f"{outpref}_counts.csv", "w") as o1:
        o1.write("Index,Token,Count\n")
        for token, count_list in token_count_list.items():
            for index, count in enumerate(count_list):
                o1.write(f"{index},{token},{count}\n")
        
    
            
    








if __name__ == "__main__":
    main()