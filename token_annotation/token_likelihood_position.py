import argparse
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
from collections import defaultdict
import math

def get_top_elements(data, N, ignore_breaks=True):
    if ignore_breaks == False:
        flat = [item for sublist in data for item in sublist if item != None]
    else:
        flat = [item for sublist in data for item in sublist if item != None and item != "break"]
    total_counts = Counter(flat)
    return total_counts, set(elem for elem, _ in total_counts.most_common(N))

def build_total_counts(total_counts, top_elements, ignore_others=False):
    top_counts = {e: total_counts[e] for e in top_elements}
    if not ignore_others:
        other_count = sum(v for k, v in total_counts.items() if k not in top_elements)
        top_counts["Other"] = other_count
    return top_counts

def plot_total_counts(total_counts, outpref, top_n, ignore_others=False):
    plt.figure(figsize=(6, 4))
    elements = list(total_counts.keys())
    values = [total_counts[e] for e in elements]

    colors = plt.cm.tab10(np.linspace(0, 1, len(elements)))
    plt.bar(elements, values, color=colors)
    plt.xlabel("Element")
    plt.ylabel("Count")
    plt.tight_layout()
    fname = f"{outpref}_total_top{top_n}" + ("_noother" if ignore_others else "") + ".png"
    plt.savefig(fname, dpi=300)
    plt.close()

def build_position_profiles_dict(data, top_elements, ignore_others=False):
    """
    Build a dictionary where each key is an element label and each value is a list of counts per position.
    """
    data_T = list(zip(*data))  # Transpose so positions are outermost
    profiles = {label: [] for label in top_elements}
    if not ignore_others:
        profiles["Other"] = []

    for pos in data_T:
        count = Counter(pos)
        for label in top_elements:
            profiles[label].append(count.get(label, 0))
        if not ignore_others:
            other_count = sum(v for k, v in count.items() if k not in top_elements)
            profiles["Other"].append(other_count)

    return profiles

def plot_profiles_from_dict(profiles, outpref, top_n, ignore_others=False, split_index=None):
    labels = list(profiles.keys())
    colors = plt.cm.tab10(np.linspace(0, 1, len(labels)))
    
    num_positions = len(next(iter(profiles.values())))
    
    # X-axis as integers, customized around split point
    if split_index is not None:
        xvals = list(range(-split_index, 0)) + list(range(0, num_positions - split_index))
    else:
        xvals = list(range(num_positions))

    # Line plot
    plt.figure(figsize=(10, 6))
    for i, label in enumerate(labels):
        plt.plot(xvals, profiles[label], label=label, color=colors[i])
    if split_index is not None:
        plt.axvline(0, color='black', linestyle='--', linewidth=1, label='Reference (0)')
    max_ticks = 20
    tick_interval = max(1, math.ceil(len(xvals) / max_ticks))
    plt.xticks(xvals[::tick_interval])
    plt.xlabel("Position")
    plt.ylabel("Count")
    plt.legend(title="Element", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    fname = f"{outpref}_profile_top{top_n}_line" + ("_noother" if ignore_others else "") + ".png"
    plt.savefig(fname, dpi=300)
    plt.close()

    # Stacked bar plot
    bottom = np.zeros(num_positions)
    plt.figure(figsize=(10, 6))
    for i, label in enumerate(labels):
        plt.bar(xvals, profiles[label], bottom=bottom, label=label, color=colors[i])
        bottom += np.array(profiles[label])
    if split_index is not None:
        plt.axvline(0, color='black', linestyle='--', linewidth=1, label='Reference (0)')
    max_ticks = 20
    tick_interval = max(1, math.ceil(len(xvals) / max_ticks))
    plt.xticks(xvals[::tick_interval])
    plt.xlabel("Position")
    plt.ylabel("Count")
    plt.legend(title="Element", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    fname = f"{outpref}_profile_top{top_n}_stacked" + ("_noother" if ignore_others else "") + ".png"
    plt.savefig(fname, dpi=300)
    plt.close()

def get_options():
    description = "Get the position profile of the highest token likelihood position."
    parser = argparse.ArgumentParser(description=description,
                                        prog='python token_likelihood_position.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--token-likelhoods',
                    required=True,
                    help='Token likehlihoods generated from compute_token_likelihood.py')
    IO.add_argument('--genomes',
                    required=True,
                    help='Input genomes used for compute_token_likelihood.py')
    IO.add_argument('--topN',
                    default=5,
                    type=int,
                    help='Number of elements to show in profile. Default = 5')
    IO.add_argument('--profile-size',
                    default=10,
                    type=int,
                    help='How many genes upstream and downstream to generate profile from. Default = 10')
    IO.add_argument('--ignore-breaks',
                    default=False,
                    action="store_true",
                    help='Ignore contig breaks in plotting.')
    IO.add_argument('--ignore-others',
                    default=False,
                    action="store_true",
                    help='Ignore genes outside of topN frequency in regions.')
    IO.add_argument('--outpref',
                default="output",
                help='Output prefix. Default = "output"')
    return parser.parse_args()

def main():
    options = get_options()
    token_likelihoods = options.token_likelhoods
    genomes = options.genomes
    profile_size = options.profile_size
    outpref = options.outpref
    topN = options.topN
    ignore_breaks = options.ignore_breaks
    ignore_others = options.ignore_others

    likehlihoods_df = pd.read_csv(token_likelihoods, index_col=0, sep='\t')
    likehlihoods_df = likehlihoods_df.drop(likehlihoods_df.columns[0], axis=1)
    likehlihoods_df.replace(0, np.nan, inplace=True)
    #print(likehlihoods_df)

    maxValueIndex = likehlihoods_df.idxmax(axis=1)
    maxValues = likehlihoods_df.max(axis=1)

    #print(maxValueIndex)
    #print(maxValues)

    file_idx = 0
    profile_list = []
    max_index = len(likehlihoods_df.columns)
    with open(genomes, "r") as f1:
        while True:
            line = f1.readline()
            if not line:
                break
            
            
            # check if name present
            if "\t" in line:
                line = line.split("\t")[-1]

            line = line.rstrip().split(" ")
            #print(line)

            # get range of genes around insertion site
            ind_post = int(maxValueIndex[file_idx])
            ind_post_range = [x if x >= 0 and x < len(line) else None for x in range(ind_post - profile_size, ind_post + profile_size)]
            #print(ind_post_range)

            pos_list = []
            for pos_idx in ind_post_range:
                if pos_idx != None:
                    #print(pos_idx)
                    pos_list.append(line[pos_idx] if line[pos_idx] != "_" else "break")
                else:
                    pos_list.append(None)

            profile_list.append(pos_list)

            file_idx += 1

    print(profile_list)

    total_counts, top_elements = get_top_elements(profile_list, topN, ignore_breaks)
    top_total_counts = build_total_counts(total_counts, top_elements, ignore_others)
    #print(total_counts)
    #print(top_elements)
    #print(top_total_counts)
    plot_total_counts(top_total_counts, outpref, topN, ignore_others)

    profiles = build_position_profiles_dict(profile_list, top_elements, ignore_others)
    #print(profiles)
    plot_profiles_from_dict(profiles, outpref, topN, ignore_others, profile_size)

if __name__ == "__main__":
    main()