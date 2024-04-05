import argparse
import sys
import os
from tokenizers import Tokenizer
import torch
from contextlib import nullcontext
import pickle
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


def get_options():
    description = "Generates genomes from LLM"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python tokenise_clusters.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--outdir',
                    default="predictions",
                    help='Output directory. Default = "predictions"')
    IO.add_argument('--reps',
                    required=True,
                    help='Output .pkl from group_sequences.py')
    IO.add_argument('--syntenyLLM',
                    required=True,
                    help='Path to output from nanoGPT trained on gene synteny.')
    IO.add_argument('--synteny_tokeniser',
                    required=True,
                    help='Path to synteny tokeniser.bin')
    IO.add_argument('--geneLLM',
                required=True,
                help='Path to output from nanoGPT trained on gene sequences.')
    IO.add_argument('--gene_tokeniser',
        required=True,
        help='Path to gene tokens.bin')
    IO.add_argument('--nanoGPT',
            required=True,
            help='Path to nanoGPT directory')
    IO.add_argument('--num_samples',
                    type=int,
                    default=1,
                    help='Number of genomes to generate. Default = 1')
    IO.add_argument('--max_genome_size',
                    type=int,
                    default=3000,
                    help='Maximum number of gene tokens to generate. Default = 3000')
    IO.add_argument('--min_genome_size',
                type=int,
                default=100,
                help='Minumum number of gene tokens to generate. Default = 100')
    IO.add_argument('--max_gene_prop',
                    type=float,
                    default=1.5,
                    help='Maximum length proportion simulated gene sequence can be of representative gene. Default = 1.5')
    IO.add_argument('--temperature',
                type=int,
                default=0.8,
                help='Randomness of sampling. 1.0 = no change, < 1.0 = less random, > 1.0 = more random. Default = 0.8')
    IO.add_argument('--top_k',
                type=int,
                default=200,
                help='Retain only the top_k most likely tokens, clamp others to have 0 probability. Default = 0.8')
    IO.add_argument('--seed',
                    type=int,
                    default=None,
                    help='Seed for sequence sampling. Default = None')
    IO.add_argument('--device',
                    default="cuda",
                    help="Device to use. Can be one of 'cpu', 'cuda', 'cuda:0', 'cuda:1', etc. Default = 'cuda'")
    IO.add_argument('--no_compile',
                    action="store_false",
                    default=True,
                    help='Do not compile model using Pytorch2 to increase speed.')
    
    

    return parser.parse_args()

def load_LLM(seed, ckpt_path, device, compile):
    from model import GPTConfig, GPT
    
    if seed != None:
        torch.manual_seed(seed)
        torch.cuda.manual_seed(seed)
    torch.backends.cuda.matmul.allow_tf32 = True # allow tf32 on matmul
    torch.backends.cudnn.allow_tf32 = True # allow tf32 on cudnn

    # model
    checkpoint = torch.load(ckpt_path, map_location=device)
    gptconf = GPTConfig(**checkpoint['model_args'])
    model = GPT(gptconf)
    state_dict = checkpoint['model']
    unwanted_prefix = '_orig_mod.'
    for k,v in list(state_dict.items()):
        if k.startswith(unwanted_prefix):
            state_dict[k[len(unwanted_prefix):]] = state_dict.pop(k)
    model.load_state_dict(state_dict)

    model.eval()
    model.to(device)
    if compile:
        model = torch.compile(model) # requires PyTorch 2.0 (optional)
    
    return model

def sample_LLM(model, device, tokenizer, max_new_tokens, temperature, top_k, start, stop, num_samples, ctx):
    # encode the beginning of the prompt
    start_ids = tokenizer.encode(start).ids
    if stop != "False":
        stop = tokenizer.encode(stop).ids[0]
    x = (torch.tensor(start_ids, dtype=torch.long, device=device)[None, ...])

    # run generation of single sequence
    generated_sequences = []
    with torch.no_grad():
        with ctx:
            for _ in range(num_samples):
                y = model.generate(x, max_new_tokens, temperature=temperature, top_k=top_k, stop = stop)
                generated_sequences.append(tokenizer.decode(y[0].tolist()))
    
    return generated_sequences

def main():
    # parse options
    #options = get_options()
    #outdir = options.outdir
    #synteny_LLM = options.synteny_LLM
    #synteny_tokeninser_path = options.synteny_tokeniser
    #gene_LLM = options.gene_LLM
    #gene_tokeniser_path = options.gene_tokeniser
    #nanoGPT_path = options.nanoGPT
    #clusters = options.clusters
    #device = options.device
    #num_samples = options.num_samples
    #temperature = options.temperature
    #top_k = options.top_k
    #seed = options.seed
    #max_genome_size = options.max_genome_size
    #max_gene_prop = options.max_gene_prop
    #compile = options.no_compile
    #min_genome_size = options.min_genome_size

    #for debugging
    outdir = "LLM_test"
    synteny_LLM = "/media/mirrored-hdd/shorsfield/jobs/pangenome_LLM/models/synteny_char/ckpt.pt"
    synteny_tokeninser_path = "/home/shorsfield/software/pangenome_LLM/data/synteny_char/tokens.bin"
    gene_LLM = "/media/mirrored-hdd/shorsfield/jobs/pangenome_LLM/models/gene_char/ckpt.pt"
    gene_tokeniser_path = "/home/shorsfield/software/pangenome_LLM/data/gene_char/tokens.bin"
    nanoGPT_path = "/home/shorsfield/software/nanoGPT"
    reps = "/home/shorsfield/software/pangenome_LLM/grouped_genes.pkl"
    device = "cuda"
    num_samples = 1
    temperature = 0.7
    top_k = 200
    seed = None
    max_genome_size = 100
    min_genome_size = 100
    max_gene_prop = 1.5
    compile = True

    # ensure min_genome_size is less than max_genome_size
    min_genome_size = min([min_genome_size, max_genome_size])

    dtype = 'bfloat16' if torch.cuda.is_available() and torch.cuda.is_bf16_supported() else 'float16' # 'float32' or 'bfloat16' or 'float16'
    device_type = 'cuda' if 'cuda' in device else 'cpu' # for later use in torch.autocast
    ptdtype = {'float32': torch.float32, 'bfloat16': torch.bfloat16, 'float16': torch.float16}[dtype]
    ctx = nullcontext() if device_type == 'cpu' else torch.amp.autocast(device_type=device_type, dtype=ptdtype)

    # import nanoGPT functions
    sys.path.append(nanoGPT_path)

    # load tokenizers and model
    special_tokens = ("[START]", "[END]")
    tokenizer = Tokenizer.from_file(synteny_tokeninser_path)
    synteny_model = load_LLM(seed, synteny_LLM, device, compile)

    # generate all genomes first, ensure all genomes are larger than min_genome_size
    print("Generating {} genomes...".format(num_samples))
    total_genomes = 0
    pred_genomes = []
    while total_genomes < num_samples:
        sim_genomes = sample_LLM(synteny_model, device, tokenizer, max_genome_size, temperature, top_k, special_tokens[0], special_tokens[1], num_samples - total_genomes, ctx)
        
        for genome in sim_genomes:
            genome = genome.replace(special_tokens[0], "").replace(special_tokens[1], "")
            genome = genome.split(" ")
            genome = [x for x in genome if x != ""]
            print(genome)
            genome_length = len(genome)
            
            print("Genome length: {}".format(genome_length))
            if genome_length >= min_genome_size:
                pred_genomes.append(genome)
        
        total_genomes += len(pred_genomes)

    
    #print(pred_genomes)
    del synteny_model

    # generate representative gene dictionary
    with (open(reps, "rb")) as f:
        reps_seq_dict = pickle.load(f)

    # sample from gene model
    gene_model = load_LLM(seed, gene_LLM, device, compile)
    tokeniser = Tokenizer.from_file(gene_tokeniser_path)
    
    # iterate through each genome, generating new gene sequence
    pred_genome_sequences = []
    genome_token_sequences = []
    print("Generating gene sequences...")
    for genome in tqdm(pred_genomes):
        genes_sequences = []
        gene_token_sequences = []
        for gene in tqdm(genome, leave=False):
            # query model
            gene_token = "_"
            if gene not in special_tokens and gene != "_":
                gene_token = int(gene)
                sequence =  reps_seq_dict[abs(gene_token)]
                sequence_length = len(sequence)
                sequence = "[START] " + sequence + " [SEP]"
                #print(sequence)

                pred_sequence = sample_LLM(gene_model, device, tokeniser, round(sequence_length * max_gene_prop), temperature, top_k, sequence, special_tokens[1], 1, ctx)
                #print(pred_sequence)

                # parse newly predicted sequences
                pred_sequence = pred_sequence[0]
                pred_sequence = "".join(pred_sequence)
                pred_sequence = pred_sequence.replace(special_tokens[1], "")
                pred_sequence = pred_sequence.replace(" ", "")
                #print(pred_sequence)
                pred_sequence = pred_sequence[sequence_length:]
                #print(pred_sequence)

                # reverse complement if required
                if gene_token < 0:
                    seq = Seq(pred_sequence)
                    pred_sequence = str(seq.reverse_complement())

                genes_sequences.append(pred_sequence)
            elif gene == "_":
                genes_sequences.append("_")
            
            gene_token_sequences.append(gene_token)
        
        # collect all new gene sequences
        pred_genome_sequences.append(genes_sequences)
        genome_token_sequences.append(gene_token_sequences)

    # generate output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # save list of genome sequences
    data = (genome_token_sequences, pred_genome_sequences)
    with open(outdir + "/predicted_sequences.pkl", "wb") as f:
        pickle.dump(data, f)  

    # concatenate gene sequences per genome
    for genome_idx, genome in enumerate(pred_genome_sequences):
        genome_seq = "".join(genome)

        # generate contigs
        genome_seq = genome_seq.split("_")
        genome_seq = [x for x in genome_seq if x != ""]

        output_sequences = []
        for contig_index, contig in enumerate(genome_seq):
            output_sequences.append(SeqRecord(Seq(contig), id="Contig_" + str(contig_index), description=""))
        
        SeqIO.write(output_sequences, outdir + "/genome_" + str(genome_idx) + ".fasta", "fasta")
    

if __name__ == "__main__":
    main()