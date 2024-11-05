import argparse
import matplotlib.pyplot as plt

def parse_args():
    """
    Parse command-line arguments.

    This function parses the command-line arguments provided by the user and returns
    a Namespace object containing the parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Token prediction with a Transformer or Reformer model.")
    parser.add_argument("--infile", type=str, required=True, help="Path to SHAP output file.")
    parser.add_argument("--tokens", type=str, required=True, help="Comma separated list of tokens to plot SHAP values for.")
    parser.add_argument("--locus", type=int, required=True, help="Position used for SHAP calculation")
    parser.add_argument("--outpref", default="plot", help="Output prefix")

    args = parser.parse_args()

    return args

def main():

    args = parse_args()
    infile = args.infile
    tokens = args.tokens.split(",")
    locus = args.locus
    outpref = args.outpref

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
                split_line = [float(x) for x in split_line[1:]]
                tokens_dict[token] = split_line
    
    # for each token, plot the change in SHAP value, ignoring locus
    for token, SHAPs_list in tokens_dict.items():
        plt.figure(figsize=(8, 6))

        # Scatter plot
        token_loci = range(len(header[:-1]))
        #plt.scatter(token_loci, SHAPs_list[:-1], color='blue')

        # Line plot connecting consecutive points, ignore base value
        plt.plot(token_loci, SHAPs_list[:-1], color='orange', linestyle='-', linewidth=1)

        plt.axvline(x=locus, color='gray', linestyle='-', linewidth=2.0, label='Reference Locus')
        plt.axhline(y=0, color='black', linestyle='--', linewidth=1.0)


        # Adding labels and legend
        plt.xlabel("Locus Token")
        plt.ylabel("SHAP value")
        plt.legend()

        plt.savefig(outpref + "_geneid_" + str(token) + "_pos_" + str(locus) + ".png", dpi=300, bbox_inches="tight")
    
if __name__ == "__main__":
    main()