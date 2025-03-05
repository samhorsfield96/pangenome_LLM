import argparse
import pandas as pd

def parse_dates(date_str):
    if len(date_str) == 4:  # Only year
        return date_str + '-01-01'
    elif len(date_str) == 7:  # Year and month
        return date_str + '-01'
    else:
        return date_str

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
    parser.add_argument("--date-type", type=str, default="collection_date", choices=['collection_date', 'last_updated'], help="Date on which to sample.")

    args = parser.parse_args()

    return args

def main():
    options = get_options()
    infile = options.infile
    outpref = options.outpref
    downsample = options.downsample
    split_by = options.split_by
    date_type = options.date_type

    df = pd.read_csv(infile, header = 0, sep = "\t")
    parsed_df = df[["sample_accession", "country", "collection_date", "last_updated", "scientific_name", "taxonomic_classification"]]
    del df
    
    if downsample != None:
        downsample_df =  pd.read_csv(downsample, header = None, sep = "\t")
        downsample_df.columns = ['sample_accession']

        parsed_df = parsed_df.merge(downsample_df, how="inner", on=["sample_accession"])

    #print(parsed_df)

    if split_by != "none":
        parsed_df[date_type] = parsed_df[date_type].apply(parse_dates)
        parsed_df[date_type] = pd.to_datetime(parsed_df[date_type], errors='coerce')
        parsed_df = parsed_df.dropna(subset=[date_type])

    #print(parsed_df[date_type])

    if split_by == "year":
        parsed_df['parsed_date'] = parsed_df[date_type].dt.year
    elif split_by == "month":
        parsed_df['parsed_date'] = parsed_df[date_type].dt.to_period('M')
    elif split_by == "day":
        parsed_df['parsed_date'] = parsed_df[date_type].dt.date

    if split_by != "none":
        for date, group in parsed_df.groupby('parsed_date'):
            group.to_csv(f'{outpref}_{date}.tsv', index=False)
    else:
        parsed_df.to_csv(f'{outpref}_all.tsv', sep='\t')


if __name__ == "__main__":
    main()