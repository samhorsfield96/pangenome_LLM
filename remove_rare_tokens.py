import argparse
import pickle
from collections import Counter

def parse_args():
    """
    Parse command-line arguments.

    This function parses the command-line arguments provided by the user and returns
    a Namespace object containing the parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Assign blanket token to genes with specific count.")
    parser.add_argument("--min", type=int, default=1, help="Minimum count of token to be included. Values equal or below will be replaced")
    parser.add_argument("--infile", required=True, help="Input file.")
    parser.add_argument("--counts", default=None, help="Path to counts generated from previous run, used for updating of existing genomes.")
    parser.add_argument("--outpref", type=str, default="genomes", help="Output prefix. Default = genomes")

    args = parser.parse_args()

    return args

def main():
    args = parse_args()
    min_count = args.min
    infile = args.infile
    tokens_path = args.tokens
    outpref = args.outpref

    gene_presence_total = Counter()
    if tokens_path != None:
        gene_presence_total = pickle.load(tokens_path)
    else:
        total_gene_count = 0
        genome_count = 0
        with open(infile, "r") as f1:
            while True:
                line = f1.readline()
                if not line:
                    break
                genome_count += 1
                gene_presence_genome = {}

                # check if name present
                if "\t" in line:
                    line = line.split("\t")[-1]

                line = line.rstrip()

                split_line = line.split(" ")
                for gene in split_line:
                    if gene != "_":
                        gene_ID = abs(int(gene))

                        # only count gene once per genome, ignore paralogs
                        gene_presence_genome[gene_ID] = 1

                        # count everything including paralogs
                        # gene_frequency_total[gene_ID] += 1
                        # total_gene_count += 1
                
                # merge dictionaries
                for key, entry in gene_presence_genome.items():
                    gene_presence_total[key] += entry
        with open(outpref + "_tokens.pkl", "wb") as f:
            pickle.dump(gene_presence_total, f)


    # get maximum gene_ID
    max_gene_ID = max(gene_presence_total.keys())

    # generate tokens that are below minimum, assigning count plus the max_gene_ID to generate new token
    below_min_dict = {gene_ID:max_gene_ID + count for gene_ID, count in gene_presence_total.items() if count <= min_count}

    with open(infile, "r") as f1, open(outpref + ".txt", "w") as f2:
        while True:
            line = f1.readline()
            if not line:
                break
            
            # check if name present
            name = None
            if "\t" in line:
                name = line.split("\t")[0]
                line = line.split("\t")[-1]

            line = line.rstrip()

            split_line = line.split(" ")
            full_line = []
            for gene in split_line:
                token_to_add = gene
                if gene != "_":
                    gene_ID = abs(int(gene))
                    if gene_ID in below_min_set:
                        token_to_add = below_min_dict[gene_ID]
                full_line.append(token_to_add)
            
            if name != None:
                f2.write(name + "\t" + " ".join(full_line) + "\n")
            else:
                f2.write(" ".join(full_line) + "\n")

if __name__ == "__main__":
    main()