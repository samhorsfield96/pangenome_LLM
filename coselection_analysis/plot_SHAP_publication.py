import argparse
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

def parse_args():
    parser = argparse.ArgumentParser(description="Plot output from compute_SHAP.py.")
    parser.add_argument("--infile", type=str, required=True, help="Path to SHAP output file.")
    parser.add_argument("--tokens", type=str, required=True, help="Comma separated list of tokens to plot SHAP values for.")
    parser.add_argument("--locus", type=int, required=True, help="Position used for SHAP calculation")
    parser.add_argument("--outpref", default="plot", help="Output prefix")
    parser.add_argument("--outlier_dev", default=3.0, type=float, help="Proportion of IQR away to identify SHAP outliers. Default = 3.0")
    parser.add_argument("--min_dist", default=100, type=int, help="Minumum distance from target token to return notable associations. Defaut = 100")

    args = parser.parse_args()

    return args

def main():

    args = parse_args()
    infile = args.infile
    tokens = args.tokens.split(",")
    locus = args.locus
    outpref = args.outpref
    outlier_dev = args.outlier_dev
    min_dist = args.min_dist

    tokens_dict = {}
    token_count = 0
    header = None
    with open(infile, "r") as f1:
        header = f1.readline().rstrip().split(",")[1:]
        while True:
            # stop when all tokens found
            if token_count == len(tokens):
                break
            line = f1.readline()
            if not line:
                break

            # pull out all matching tokens
            split_line = line.rstrip().split(',')
            token = split_line[0]
            if token in tokens:
                converted_line = []
                for val in split_line[1:]:
                    try:
                        converted_line.append(float(val))
                    except (ValueError, TypeError):
                        converted_line.append(0.0)
                tokens_dict[token] = converted_line
                token_count += 1
    
    min_dist_upper = locus + min_dist
    min_dist_lower = locus - min_dist
    token_list = header[:-1]

    # contig break locations
    contig_breaks = [x for x in range(len(token_list)) if token_list[x] == "_"]
    dists_to_contig_break = [x - locus for x in contig_breaks]
    upstream_contig_break = max([x + locus for x in dists_to_contig_break if x < 0], default=0)
    downstream_contig_break = min([x + locus for x in dists_to_contig_break if x > 0], default=len(token_list))

    print(f"min_dist_lower: {min_dist_lower}")
    print(f"min_dist_upper: {min_dist_upper}")
    print(f"contig_breaks: {contig_breaks}")
    print(f"upstream_contig_break: {upstream_contig_break}")
    print(f"downstream_contig_break: {downstream_contig_break}")

    # check which is closer, a contig break of the min distance cutoff
    min_dist_lower = max(min_dist_lower, upstream_contig_break)
    min_dist_upper = min(min_dist_upper, downstream_contig_break)

    print(f"min_dist_lower: {min_dist_lower}")
    print(f"min_dist_upper: {min_dist_upper}")

    # for each token, plot the change in SHAP value, ignoring locus
    for token, SHAPs_list in tokens_dict.items():
        plt.figure(figsize=(8, 6))

        values = np.array(SHAPs_list[:-1])

        # tukey's method for outlier detection
        q1 = np.percentile(values, 25)
        q3 = np.percentile(values, 75)
        iqr = q3 - q1
        
        lower_bound = q1 - (outlier_dev * iqr)
        upper_bound = q3 + (outlier_dev * iqr)

        #outliers = values[z_scores > outlier_std]
        outlier_tokens = [(token_list[x], values[x], x) for x in range(len(token_list)) if (values[x] < lower_bound or values[x] > upper_bound) and (x < min_dist_lower or x > min_dist_upper)]

        #print(token)
        #print(outlier_tokens)

        # Scatter plot
        token_loci = range(len(token_list))
        #plt.scatter(token_loci, SHAPs_list[:-1], color='blue')

        # x and y arrays
        x = np.array(token_loci)
        y = np.array(values)

        # add contig breaks
        for contig_break in contig_breaks:
            plt.axvline(x=contig_break, color='#C0C0C0', linestyle='--', linewidth=0.5)

        # Identify outlier indices
        outlier_mask = (y > upper_bound) | (y < lower_bound)

        # Plot ONLY the outliers, no outline
        plt.scatter(
            x[outlier_mask],
            y[outlier_mask],
            color='#E64B35FF',
            s=15,
            marker='o',
            linewidths=0,     # <- no outline
            edgecolors='none', # <- no outline
            zorder=3
        )

        # Always plot the main line in grey
        plt.plot(x, y, color='black', linestyle='-', linewidth=0.5)

        plt.axvline(x=locus, color='#4DBBD5FF', linestyle='-', linewidth=1.0)
        plt.axhline(y=0, color='black', linestyle='solid', linewidth=1.0)

        #if min_dist_upper <= len(token_list):
        #    plt.axvline(x=min_dist_upper, color='red', linestyle='--')
        
        #if min_dist_lower >= 0:
        #   plt.axvline(x=min_dist_lower, color='red', linestyle='--')

        plt.axhline(y=lower_bound, color='black', linestyle='--', linewidth=1.0)
        plt.axhline(y=upper_bound, color='black', linestyle='--', linewidth=1.0)

        # Adding labels and legend
        plt.xlabel("Gene locus", fontweight='bold')
        plt.ylabel("SHAP value (log pseudolikelihood)", fontweight='bold')
        plt.legend()

        plt.savefig(outpref + "_geneid_" + str(token) + "_pos_" + str(locus) + ".png", dpi=300, bbox_inches="tight")
        plt.savefig(outpref + "_geneid_" + str(token) + "_pos_" + str(locus) + ".svg", dpi=300, bbox_inches="tight")

        with open(outpref + "_geneid_" + str(token) + "_pos_" + str(locus) + "_hits.csv", "w") as f:
           f.write("Token,SHAP,Locus\n")
           for entry in outlier_tokens:
               f.write(f"{entry[0]},{entry[1]},{entry[2]}\n")
                
    
if __name__ == "__main__":
    main()