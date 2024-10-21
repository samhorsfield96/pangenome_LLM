import argparse
import re
import networkx as nx
from collections import Counter

def get_options():
    description = "Compares synteny between simulated and generated genomes"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python synteny_accuracy.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--outpref',
                    default="comparisons",
                    help='Output prefix. Default = "comparisons"')
    IO.add_argument('--infile',
                    required=True,
                    help='Path to tokenised genome data .txt file.')
    return parser.parse_args()

def main():
    options = get_options()
    infile = options.infile
    outpref = options.outpref

    total_gene_count = 0
    gene_presence_total = Counter()
    gene_frequency_total = Counter()
    with open(infile, "r") as f:
        while True:
            line = f1.readline()
            if not line:
                break
            gene_presence_genome = {}

            line = line.rstrip()

            split_line = line.split(" ")
            for gene in split_line:
                if gene != "_":
                    gene_ID = abs(int(gene))

                    # only count gene once per genome, ignore paralogs
                    gene_presence_genome[gene_ID] = 1

                    # count everything including paralogs
                    gene_frequency_total[gene_ID] += 1
                    total_gene_count += 1
            
            # merge dictionaries
            for key, entry in gene_presence_genome.items():
                gene_presence_total[key] += entry
    
    # print to file
    total_genes_presence = len(gene_presence_total)
    with (outpref + ".txt", "w") as o:
        o.write("Gene_ID\tGenome_count\tTotal_count\tGenome_freq\tTotal_freq\n")
        for gene_ID, count in gene_presence_total.items():
            o.write("{}\t{}\t{}\t{}\t{}\n".format(gene_ID, count, gene_frequency_total[gene_ID], count/total_genes_presence, gene_frequency_total[gene_ID]/total_gene_count))


if __name__ == "__main__":
    main()