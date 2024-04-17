from pathlib import Path
import argparse
import os
import pickle
from Bio import Align
from Bio.Seq import Seq
import numpy as np
import pandas as pd

def get_options():
    description = "Generates genomes from LLM"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python tokenise_clusters.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--outpref',
                    default="comparisons",
                    help='Output prefix. Default = "comparisons"')
    IO.add_argument('--reps',
                    required=True,
                    help='Output .pkl from group_sequences.py')
    IO.add_argument('--LLM_output',
                    required=True,
                    help='Path to output .pkl from generate_genomes.py.')
    IO.add_argument('--token_data',
                    required=True,
                    help='Path to output .pkl from tokenise_clusters.py')
    return parser.parse_args()

def main():
    # options = get_options()
    # token_data = options.token_data
    # reps = options.reps
    # LLM_output = options.LLM_output
    # outpref = options.outpref
    
    token_data = "/home/shorsfield/software/pangenome_LLM/tokenised_genomes.pkl"
    reps = "/home/shorsfield/software/pangenome_LLM/grouped_genes.pkl"
    LLM_output = "/home/shorsfield/software/pangenome_LLM/LLM_out_temp_1.5/predicted_sequences.pkl"
    outpref = "temp_1.5_comparison"
    
    # load in all data
    with (open(token_data, "rb")) as f:
        gene_tokens, reps_dict = pickle.load(f)
    
    with (open(reps, "rb")) as f:
        reps_seq_dict = pickle.load(f)
    
    with (open(LLM_output, "rb")) as f:
        genome_token_sequences, pred_genome_sequences = pickle.load(f)
    
    # remove all underscores
    genome_token_sequences_iter = [[item for item in sublist if item != "_"] for sublist in genome_token_sequences]
    pred_genome_sequences_iter = [[item for item in sublist if item != "_"] for sublist in pred_genome_sequences]

    # calculate ANI 
    id_list = []
    rep_len_list = []
    seq_len_list = []
    token_list = []
    genome_id_list = []
    aligner = Align.PairwiseAligner()
    genome_id = 0
    for tokenised_genome, seq_genome in zip(genome_token_sequences_iter, pred_genome_sequences_iter):
        for gene_token, gene_seq in zip(tokenised_genome, seq_genome):
            gene_token = int(gene_token)
            rep_sequence =  reps_seq_dict[abs(gene_token)]

            rep_len_list.append(len(rep_sequence))
            seq_len_list.append(len(gene_seq))
            token_list.append(gene_token)
            genome_id_list.append(genome_id)

            # catch cases where no sequence generated
            if len(gene_seq) == 0:
                id_list.append(0.0)
                continue

            gene_seq = Seq(gene_seq)

            if gene_token < 0:
                gene_seq = gene_seq.reverse_complement()

            # take top alignment
            alignment = next(aligner.align(Seq(rep_sequence), gene_seq))

            #print("Score = %.1f" % alignment.score)
            #print(alignment)
            
            seq1, seq2 = alignment
            align1_arr = np.array(list(str(seq1)))
            align2_arr = np.array(list(str(seq2)))
            num_match = np.count_nonzero(align1_arr == align2_arr)
            id_list.append(num_match / align1_arr.size)
        genome_id += 1

    df = pd.DataFrame({'token': token_list, 'prop_id': id_list, 'rep_length': rep_len_list, 'sim_len': seq_len_list})
    
    df.to_csv(outpref + '.csv', index=False) 
    


if __name__ == "__main__":
    main()