import argparse
import re
import networkx as nx

def get_options():
    description = "Compares synteny between simulated and generated genomes"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python synteny_accuracy.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--outpref',
                    default="comparisons",
                    help='Output prefix. Default = "comparisons"')
    IO.add_argument('--real',
                    required=True,
                    help='Path to real data .txt file.')
    IO.add_argument('--sim',
                    required=True,
                    help='Path to simulated data .txt file.')
    return parser.parse_args()

def clean_genome(input_string):
    # more than one hyphen
    pattern = r'-{2,}'
    input_string = re.sub(pattern, '-', input_string)

    # hyphens on their own and leading
    pattern = r'^\s*-\s*'
    input_string = re.sub(pattern, '-', input_string)

    # hyphens on their own
    input_string = re.sub(r'\s-\s', ' ', input_string)

    # trailing hyphens
    input_string = re.sub(r'\s-\s*$', '', input_string)

    # address contig break issues
    input_string = re.sub(r'-_', '_', input_string)

    input_string = input_string.strip()

    return input_string


def main():
    #options = get_options()
    #real_genomes = options.real
    #sim_genomes = options.sim
    #outpref = options.outpref

    real_genomes = "/home/shorsfield/software/panGPT/test_prompt.txt"
    sim_genomes = "/home/shorsfield/software/panGPT/simulations_temp_1.0_BPE_tokeniser.txt"
    outpref = "sim_test.txt"

    real_genome_list = []
    with open(real_genomes, "r") as f1:
        while True:
            line = f1.readline()
            if not line:
                break
            line = line.rstrip()

            #print("Pre cleaning:")
            #print(line)
            line = clean_genome(line)

            #print("Post cleaning:")
            #print(line)

            real_genome_list.append(line.split())

    sim_genome_list = []
    with open(sim_genomes, "r") as f1:
        while True:
            line = f1.readline()
            if not line:
                break
            line = line.rstrip()

            print("Pre cleaning:")
            #print(line)
            line = clean_genome(line)

            print("Post cleaning:")
            #print(line)

            sim_genome_list.append(line.split())

    DG=nx.DiGraph()
    for index, genome in enumerate(real_genomes):
        for i in range(len(genome) - 1):
            u, v = genome[i], genome [i + 1]

            # check if nodes are present
            if u != "_":
                if not DG.has_node(abs(u)) and :
                    DG.add_node(abs(u),
                                members_real=[],
                                members_sim=[])
            


            add_cluster_u = True if not DG.has_node(u) else False
            add_cluster_v = True if not DG.has_node(u) else False

            # encode directionality, if negative reverse orientation
            
            DG.add_node(abs(v), label='real')
            DG.add_edge(u, v, label='sim')



if __name__ == "__main__":
    main()