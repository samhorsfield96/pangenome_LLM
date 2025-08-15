import argparse
from pathlib import Path

def parse_args():
    """
    Parse command-line arguments.

    This function parses the command-line arguments provided by the user and returns
    a Namespace object containing the parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Pull out filenames that match between two files.")
    parser.add_argument("--index", type=str, required=True, help="Path to the index file")
    parser.add_argument("--target", type=str, required=True, help="Path to the target file to search in")
    parser.add_argument("--outfile", type=str, default="output.txt", help="Output filename. Default = output.txt")

    args = parser.parse_args()

    return args

def main():
    args = parse_args()
    index = args.index
    target = args.target
    outfile = args.outfile

    index_set = set()
    with open(index, "r") as f1:
        for line in f1:
            parsed_line = line.rstrip()
            # remove separator stop if added
            if parsed_line[-1] == ".":
                parsed_line = parsed_line[:-1]
            index_set.add(parsed_line)
    
    with open(outfile, "w") as o1:
        with open(target, 'r') as f2:
            for line in f2:
                parsed_line = line.rstrip()
                name = Path(parsed_line).name.split('.')[0]

                if name in index_set:
                    o1.write(line)




if __name__ == "__main__":
    main()