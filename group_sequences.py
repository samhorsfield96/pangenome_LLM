from pathlib import Path
import argparse
import os
import pickle
from Bio import SeqIO

def get_options():
    description = "Groups clusters into single fasta files and assigns token."
    parser = argparse.ArgumentParser(description=description,
                                        prog='python group_sequences.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--pickle',
                    required=True,
                    help='Pickle file from tokenise_clusters')
    IO.add_argument('--fasta',
                    required=True,
                    help='Nucleotide fasta file containing all gene sequences concetenated into single file.')
    IO.add_argument('--outpref',
                default="grouped_genes",
                help='Output prefix.')

    return parser.parse_args()

def get_gff(directory):
    files = (str(p.resolve()) for p in Path(directory).rglob("*.gff3"))
    yield from files

def main():
    options = get_options()
    pkl = options.pickle
    fasta = options.fasta
    outpref = options.outpred

    #pkl = "/home/shorsfield/software/pangenome_LLM/tokenised_genomes.pkl"
    #fasta = "/media/mirrored-hdd/shorsfield/jobs/pangenome_LLM/all_seqs.ffn"
    #outpref = "grouped_genes"

    with (open(pkl, "rb")) as f:
        gene_tokens, reps_dict = pickle.load(f)

    
    reps_seq_dict = {}
    gene_dict = {}
    fasta_sequences = SeqIO.parse(open(fasta),'fasta')
    for fasta in fasta_sequences:
        gene_id, sequence = fasta.id, str(fasta.seq)

        gene_token = gene_tokens.get(gene_id, None)
        if gene_token != None:
            curr_rep = reps_dict[gene_token]
            if gene_id == curr_rep:
                reps_seq_dict[gene_token] = sequence
            else:
                if gene_token not in gene_dict:
                    gene_dict[gene_token] = []
                gene_dict[gene_token].append(sequence)
    

    with open(outpref + ".txt", "w") as o:
        for gene_token, gene_list in gene_dict.items():
            # write reprsentative, comma and each gene
            rep_seq = reps_seq_dict[gene_token]
            for gene in gene_list:
                o.write(rep_seq + "," + gene + "\n")                 

if __name__ == "__main__":
    main()