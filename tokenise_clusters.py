from pathlib import Path
import argparse
import os
import pickle

def get_options():
    description = "Tokenises gene clusters"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python tokenise_clusters.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--tool',
                required=True,
                help='Tool used for clustering. One of mmseqs2 or CD-Hit')
    IO.add_argument('--gffs',
                    required=True,
                    help='Input directory containing gff files. Will search recursively')
    IO.add_argument('--clusters',
                    required=True,
                    help='Cluster file.')
    IO.add_argument('--outpref',
                default="tokenised_genomes",
                help='Output prefix.')

    return parser.parse_args()

def get_gff(directory):
    files = (str(p.resolve()) for p in Path(directory).rglob("*.gff3"))
    yield from files

def main():
    options = get_options()
    tool = options.tool
    gff_dir = options.gffs
    cluster_file = options.clusters
    outpref = options.outpref

    #tool = "mmseqs2"
    #gff_dir = "/media/mirrored-hdd/shorsfield/jobs/pangenome_LLM/bakta"
    #outpref = "tokenised_genomes"
    #cluster_file = "/media/mirrored-hdd/shorsfield/jobs/pangenome_LLM/mmseqs_id60_len60_cluster_sorted.tsv"

    # dictionary of representative sequences and their token
    reps_dict = {}

    # dictionary mapping each gene to a given cluster token
    gene_tokens = {}

    token = -1
    if tool == "mmseqs2":
        current_rep = None
        with open(cluster_file, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break

                split_line = line.rstrip().split("\t")
                rep = split_line[0]
                seq = split_line[1]

                # new cluster, increment token
                if current_rep != rep:
                    current_rep = rep
                    token += 1
                    reps_dict[token] = current_rep
                
                # add sequence to cluster
                gene_tokens[seq] = token

    #print(len(reps_dict))

    # generate list of integers to represent genome
    genome_list = []
    for gff in get_gff(gff_dir):
        basename = os.path.basename(gff)
        with open(gff, "r") as f:
            tokenised_genome = []
            while True:
                line = f.readline()
                if not line or line == "##FASTA\n":
                    break
                
                # skip commented lines
                if line[0] == "#":
                    continue

                split_line = line.rstrip().split("\t")
                type = split_line[2]
                if type == "region":
                    # add space between contigs as synteny is unknown
                    if len(tokenised_genome) > 0:
                        tokenised_genome.append("_")
                else:
                    gene_strand = True if split_line[6] == "+" else False
                    gene_id = split_line[-1].split(";")[0].replace("ID=", "")
                    gene_token = gene_tokens.get(gene_id, None)

                    if gene_token != None:
                        # multiply by minus 1 for negative strand
                        if not gene_strand:
                            gene_token *= -1
                        
                        tokenised_genome.append(str(gene_token))
        
        tokenised_genome_str = " ".join(tokenised_genome)
        genome_list.append((basename, tokenised_genome_str))
    
    with open(outpref + ".txt", "w") as o:
        for genome_name, tokenised_genome_str in genome_list:
            o.write(genome_name + "\t" + tokenised_genome_str + "\n")
    
    # save data as pickle
    data = (gene_tokens, reps_dict)

    with open(outpref + ".pkl", "wb") as f:
        pickle.dump(data, f)                  

if __name__ == "__main__":
    main()