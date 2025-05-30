from pathlib import Path
import argparse
import os
import pickle
#from multiprocessing import Pool, Manager
#from functools import partial
import rocksdb
import math

def get_options():
    description = "Tokenises gene clusters"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python tokenise_clusters.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--gffs',
                    default=None,
                    help='File containing gff files to search for, one absolute path per line. [Default = None]')
    IO.add_argument('--clusters',
                    required=True,
                    help='Cluster file.')
    IO.add_argument('--outpref',
                default="tokenised_genomes",
                help='Output prefix.')
    IO.add_argument('--db',
            default=None,
            help='Path to previous gene tokens db. [Default=None]')   
    IO.add_argument('--process_id',
        default=1,
        type=int,
        help='Process ID. [Default = 1]')

    return parser.parse_args()

def get_gff(directory):
    files = (str(p.resolve()) for p in Path(directory).rglob("*.gff*"))
    return list(files)

def chunks(l, n):
    """Yield n number of striped chunks from l."""
    for i in range(0, n):
        yield l[i::n]

def generate_gene_id(unsplit_id):
    name = unsplit_id.split("SAM")[1].split("_")[0].split(".")[0]
    split_gene = unsplit_id.split("_")
    gene_id = split_gene[-1]
    contig_id = split_gene[-2][-5:]

    gene_name = name + "_" + contig_id + "_" + gene_id

    return gene_name

def tokenise_gff(index, gff_list, outpref, gene_tokens):
    with open(outpref + "_batch_" + str(index) + ".txt", "w") as o:
        for gff in gff_list:
            basename = os.path.basename(gff)
            with open(gff, "r") as f:
                tokenised_genome = []
                current_contig = None
                while True:
                    line = f.readline()
                    if not line or line == "##FASTA\n":
                        break
                    
                    # skip commented lines
                    if line[0] == "#":
                        continue

                    split_line = line.rstrip().split("\t")
                    #type = split_line[2]
                    # if type == "region":
                    #     # add space between contigs as synteny is unknown
                    #     if len(tokenised_genome) > 0:
                    #         tokenised_genome.append("_")
                    # else:
                    
                    gene_strand = True if split_line[6] == "+" else False
                    split_gene_id = split_line[-1].split(";")[0].replace("ID=", "").split("_")

                    # get contig name 
                    name = split_line[0].replace("SAM", "").split(".contig") #basename.split("SAM")[1].split("_")[0].split(".")[0]
                    
                    # deal with issue splitting at contig
                    try:
                        contig_ID = name[1]
                    except:
                        print(f"Can't split {split_line[0]}")
                        continue
                    
                    gene_ID = split_gene_id[1]

                    # add contig end
                    if contig_ID != current_contig:
                        if len(tokenised_genome) > 0:
                            tokenised_genome.append("_")
                        current_contig = contig_ID

                    # build gene id to search in dictionary
                    gene_name = name[0] + "_" + contig_ID + "_" + gene_ID
                    #print(gene_name)

                    gene_token = gene_tokens.get(gene_name.encode())
                    if gene_token is not None:
                        gene_token = gene_token.decode()
                        # multiply by minus 1 for negative strand
                        if not gene_strand:
                            gene_token = "-" + gene_token
                        
                        tokenised_genome.append(str(gene_token))
                    else:
                        print(f"{gene_name} not found")
            
            tokenised_genome_str = " ".join(tokenised_genome)
            o.write(basename + "\t" + tokenised_genome_str + "\n")
    return index

def main():
    options = get_options()
    gff_dir = options.gffs
    cluster_file = options.clusters
    outpref = options.outpref
    gene_tokens_db = options.db
    process_id = options.process_id

    #gff_dir = "/media/mirrored-hdd/shorsfield/jobs/pangenome_LLM/bakta"
    #outpref = "tokenised_genomes"
    #cluster_file = "/media/mirrored-hdd/shorsfield/jobs/pangenome_LLM/mmseqs_id60_len60_cluster_sorted.tsv"

    opts = rocksdb.Options()
    opts.max_open_files = 300000000
    opts.max_bytes_for_level_base = 209715200 #(default = 10485760)
    opts.target_file_size_base = math.ceil(opts.max_bytes_for_level_base / 10) #(default = 2097152)
    opts.target_file_size_multiplier = 2 #(default = 1)

    # dictionary of representative sequences and their token
    if gene_tokens_db is None:
        reps_dict = {}
        rep_to_token = {}

        # dictionary mapping each gene to a given cluster token
        gene_tokens_db = outpref + "_gene_tokens.db"
        opts.create_if_missing = True
        gene_tokens = rocksdb.DB(gene_tokens_db, opts)

        #start at 0 as cannot assign negative 0
        token = 0
        current_rep = None
        print("Generating token dictionaries...")
        counter = 0
        with open(cluster_file, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break

                split_line = line.rstrip().split("\t")
                split_rep = split_line[0]
                split_seq = split_line[1]

                rep = generate_gene_id(split_rep)
                seq = generate_gene_id(split_seq)
                #print(rep)
                #print(seq)

                # allows use of non-sorted list
                if rep not in rep_to_token:
                    token += 1
                    rep_to_token[rep] = token
                
                current_token = rep_to_token[rep]
                reps_dict[current_token] = rep
                
                # add sequence to cluster
                gene_tokens.put(seq.encode(), str(current_token).encode())
                counter += 1
                if counter % 10000000 == 0:
                    print("At index: {}".format(counter))

        # save data as pickle
        print("Saving token dictionaries...")

        with open(outpref + "_reps.pkl", "wb") as f:
            pickle.dump(reps_dict, f)
        
        del reps_dict
    
        print("Saved token dictionaries.")
    else:
        print("Loading db...")
        gene_tokens = rocksdb.DB(gene_tokens_db, opts, read_only=True)
    
        print("Generating tokenised genomes...")

        # generate list of integers to represent genome
        files_list = []
        #files_list = get_gff(gff_dir)
        #file_list_chunks = chunks(files_list, num_processes)

        with open(gff_dir, "r") as o:
            while True:
                line = o.readline()
                if not line:
                    break
                files_list.append(line.rstrip())

        index = tokenise_gff(process_id, files_list, outpref=outpref, gene_tokens=gene_tokens)
        print("Finished batch: {}".format(str(process_id)))

if __name__ == "__main__":
    main()
