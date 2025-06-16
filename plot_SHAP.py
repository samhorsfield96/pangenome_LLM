import argparse
import matplotlib.pyplot as plt
import numpy as np

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
    header = None
    with open(infile, "r") as f1:
        header = f1.readline().rstrip().split(",")[1:]
        while True:
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
    
    min_dist_upper = locus + min_dist
    min_dist_lower = locus - min_dist

    # for each token, plot the change in SHAP value, ignoring locus
    for token, SHAPs_list in tokens_dict.items():
        plt.figure(figsize=(8, 6))

        values = np.array(SHAPs_list[:-1])
        token_list = header[:-1]

        # tukey's method for outlier detection
        q1 = np.percentile(values, 25)
        q3 = np.percentile(values, 75)
        iqr = q3 - q1
        
        lower_bound = q1 - (outlier_dev * iqr)
        upper_bound = q3 + (outlier_dev * iqr)

        #outliers = values[z_scores > outlier_std]
        outlier_tokens = [(token_list[x], values[x], x) for x in range(len(token_list)) if (values[x] < lower_bound or values[x] > upper_bound) and abs(x - locus) >= min_dist]

        #print(token)
        #print(outlier_tokens)

        # Scatter plot
        token_loci = range(len(token_list))
        #plt.scatter(token_loci, SHAPs_list[:-1], color='blue')

        # Line plot connecting consecutive points, ignore base value
        plt.plot(token_loci, values, color='orange', linestyle='-', linewidth=1)

        plt.axvline(x=locus, color='gray', linestyle='-', linewidth=2.0, label='Reference Locus')
        plt.axhline(y=0, color='black', linestyle='--', linewidth=1.0)

        if min_dist_upper <= len(token_list):
            plt.axvline(x=min_dist_upper, color='red', linestyle='--')
        
        if min_dist_lower >= 0:
            plt.axvline(x=min_dist_lower, color='red', linestyle='--')

        plt.axhline(y=lower_bound, color='red', linestyle='--')
        plt.axhline(y=upper_bound, color='red', linestyle='--')


        # Adding labels and legend
        plt.xlabel("Gene locus")
        plt.ylabel("SHAP value (log pseudolikelihood)")
        plt.legend()

        plt.savefig(outpref + "_geneid_" + str(token) + "_pos_" + str(locus) + ".png", dpi=300, bbox_inches="tight")

        with open(outpref + "_geneid_" + str(token) + "_pos_" + str(locus) + "_hits.csv", "w") as f:
            f.write("Token,SHAP,Locus\n")
            for entry in outlier_tokens:
                f.write(f"{entry[0]},{entry[1]},{entry[2]}\n")
                
    
if __name__ == "__main__":
    main()