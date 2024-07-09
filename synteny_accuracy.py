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

    # hyphens on their own, remove with following space
    input_string = re.sub(r'-\s', '', input_string)

    # internal hyphens, add as space
    input_string = re.sub(r'(\d)-(\d)', ' -', input_string)

    # trailing hyphens
    input_string = re.sub(r'\s-\s*$', '', input_string)

    # address contig break issues
    input_string = re.sub(r'-_', '_', input_string)

    input_string = input_string.strip()

    return input_string

def add_node(G, genome_idx, node_id, real=True):

    if node_id != "_":
        node_id = int(node_id)
        if not G.has_node(abs(node_id)):
            G.add_node(abs(node_id),
                        members_real=[],
                        members_sim=[])
    
        if real == True:
            G.nodes[abs(node_id)]['members_real'].append(genome_idx)
        else:
            G.nodes[abs(node_id)]['members_sim'].append(genome_idx)

def add_edge(G, genome_idx, u, v, real=True):

    if u != "_" and v != "_":
        u = int(u)
        v = int(v)
        u_strand = 1 if u >= 0 else 0
        v_strand = 1 if v >= 0 else 0

        # add in absolute order
        first = min(abs(u), abs(v))
        second = max(abs(u), abs(v))


        same_strand = True if u_strand == v_strand else False
        
        if not G.has_edge(first, second):
            G.add_edge(first, second, 
                            members_real=[],
                            members_sim=[],
                            strand_real=[],
                            strand_sim=[])
        
        if real == True:
            G.edges[first, second]['members_real'].append(genome_idx)
            G.edges[first, second]['strand_real'].append(same_strand)
        else:
            G.edges[first, second]['members_sim'].append(genome_idx)
            G.edges[first, second]['strand_sim'].append(same_strand)

def main():
    #options = get_options()
    #real_genomes = options.real
    #sim_genomes = options.sim
    #outpref = options.outpref

    real_genomes = "/home/shorsfield/software/panGPT/test_prompt.txt"
    sim_genomes = "/home/shorsfield/software/panGPT/simulations_temp_1.0_BPE_tokeniser.txt"
    outpref = "sim_test"

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

            #print("Pre cleaning:")
            #print(line)
            line = clean_genome(line)

            #print("Post cleaning:")
            #print(line)

            sim_genome_list.append(line.split())

    G=nx.Graph()
    real = True
    for index, genome in enumerate(real_genome_list):
        for i in range(len(genome) - 1):
            u, v = genome[i], genome [i + 1]

            # add nodes
            add_node(G, index, u, real=real)
            add_node(G, index, v, real=real)

            # add edges, encoding directionality
            add_edge(G, index, u, v, real=real)
    
    real = False
    for index, genome in enumerate(sim_genome_list):
        print(genome)
        for i in range(len(genome) - 1):
            u, v = genome[i], genome [i + 1]

            # add nodes
            add_node(G, index, u, real=real)
            add_node(G, index, v, real=real)

            # add edges, encoding directionality
            add_edge(G, index, u, v, real=real)

    nx.write_gml(G, outpref + '.gml')

if __name__ == "__main__":
    main()