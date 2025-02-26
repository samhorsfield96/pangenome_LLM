import argparse
import pandas as pd

def get_options():
    """
    Parse command-line arguments.

    This function parses the command-line arguments provided by the user and returns
    a Namespace object containing the parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Parses samples based on their dates.")
    parser.add_argument("--infile", type=str, required=True, help="ENA metadata file.")
    parser.add_argument("--outpref", type=str, default="output", help="Output prefix. Default = 'output'")
    parser.add_argument("--downsample", type=str, default=None, help="File of single column of genome IDs to downsample by, does not have a header.")
    parser.add_argument("--split-by", type=str, default="none", choices=['year', 'month', 'day', 'none'], help="How to stratify samples by date.")

    args = parser.parse_args()

    return args

def main():
    options = get_options()
    infile = options.infile
    outpref = options.outpref
    downsample = options.downsample
    split_by = options.split_by

    df = pd.read_csv(infile, header = 0, sep = "\t")
    parsed_df = df[["sample_accession", "country", "collection_date", "scientific_name", "taxonomic_classification"]]
    del df
    
    parsed_df['collection_date'] = pd.to_datetime(parsed_df['collection_date'])

    if downsample != None:
        downsample_df =  pd.read_csv(downsample, header = False, sep = "\t")
        downsample_df.columns = ['sample_accession']

        parsed_df = parsed_df.merge(downsample_df, how="inner", on=["sample_accession"])

    if split_by == "year":
        parsed_df['parsed_date'] = parsed_df['date'].dt.year
    elif split_by == "month":
        parsed_df['parsed_date'] = parsed_df['date'].dt.to_period('M')
    elif split_by == "day":
        parsed_df['parsed_date'] = parsed_df['date'].dt.date

    if split_by != "none":
        for date, group in parsed_df.groupby('parsed_date'):
            group.to_csv(f'{outpref}_{date}.tsv', index=False)
    else:
        parsed_df.to_csv(f'{outpref}_all.tsv', sep='\t')


if __name__ == "__main__":
    main()