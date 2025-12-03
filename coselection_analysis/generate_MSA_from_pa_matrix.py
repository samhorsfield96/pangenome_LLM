import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Generate mock MSA from gene p/a matrix.")
    parser.add_argument("--infile", type=str, required=True, help="Path to gene presence/absence matrix.")
    parser.add_argument("--outpref", default="plot", help="Output prefix")

    args = parser.parse_args()

    return args

def wrap(seq, width):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def main():
    options = parse_args()
    infile = options.infile
    outpref = options.outpref

    fasta_wrap = 60

    # df: your presence/absence matrix (rows = genes, columns = genomes)


    # determine if reversed
    reverse = False
    with open(infile, "r") as input:
        header = input.readline().split(",")
        if header[0] == "genome":
            reverse = True
    
    # reverse is more computational efficient, although input is unconventional
    if reverse:
        with open(infile, "r") as input, open(outpref + ".fasta", "w") as out:
            header = input.readline()

            for line in input:
                split_line = line.rstrip().split(",")
                genome = split_line[0]
                seq = ''.join('c' if val == "1" else 'a' for val in split_line[1:])
                out.write(f">{genome}\n{wrap(seq, fasta_wrap)}\n")

    else:
        df = pd.read_csv(infile, sep=",")
        df = df.set_index(df.columns[0])

        with open(outpref + ".fasta", "w") as out:
            for genome in df.columns:
                seq = ''.join('c' if val == 1 else 'a' for val in df[genome])
                out.write(f">{genome}\n{wrap(seq, fasta_wrap)}\n")

        df = pd.read_csv(infile, sep=",")
        df = df.set_index(df.columns[0])

        for genome in df.columns:
            seq = ''.join('c' if val == 1 else 'a' for val in df[genome])
            out.write(f">{genome}\n{wrap(seq, fasta_wrap)}\n")

if __name__ == "__main__":
    main()
