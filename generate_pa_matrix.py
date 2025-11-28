import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Generate gene pa matrix from tokenised genomes.py.")
    parser.add_argument("--infile", type=str, required=True, help="Path to tokenised genome file.")
    parser.add_argument("--outpref", default="plot", help="Output prefix")

    args = parser.parse_args()

    return args

def main():
    options = parse_args()
    infile = options.infile
    outpref = options.outpref

    # Read file
    genes = {}

    with open(infile, "r") as f:
        for line in f:
            genome, vals = line.strip().split("\t")
            vals = vals.split()

            vals = [str(abs(int(v))) for v in vals if v != "_"]
            
            genes[genome] = set(vals)

    # Collect all unique genes
    all_genes = sorted({g for s in genes.values() for g in s}, key=int)

    # Build presence/absence matrix
    df = pd.DataFrame(index=all_genes, columns=genes.keys())

    for genome, gset in genes.items():
        df[genome] = df.index.map(lambda x: 1 if x in gset else 0)

    df.index.name = "gene"

    df.to_csv(outpref + ".csv")

if __name__ == "__main__":
    main()
