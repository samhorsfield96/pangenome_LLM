from collections import defaultdict
import pickle
import argparse
import sys
from tqdm import tqdm

def get_options():

    parser = argparse.ArgumentParser(description='Using output from WTBcluster, calculate pangenome statistics')

    # input options
    parser.add_argument('--query', help='Smaller of cluster files.',
                                    required=True)
    parser.add_argument('--database', help='Large of cluster files.',
                                     required=True)
    parser.add_argument('--min-cluster-size', help='Minimum size of cluster to compare',
                                     default=0)
    parser.add_argument('--outpref', help='Output prefix',
                                    required = True)

    return parser.parse_args()

def main():
    args = get_options()

    query_file = args.query
    database = args.database
    outpref = args.outpref
    min_cluster_size = args.min_cluster_size

    query_set_dict = defaultdict(set)
    query_set = set()
    # get all queries
    print("Reading queries...", file=sys.stderr)
    with open(query_file, "r") as f1:
        for line in f1:
            split_line = line.rstrip().split("\t")
            query_id = split_line[1]
            rep_id = split_line[0]
            query_set.add(query_id)

            query_set_dict[rep_id].add(query_id)
    
    # remove clusters below min cluster size
    print("Filtering clusters on size...", file=sys.stderr)
    to_remove = set()
    if min_cluster_size > 0:
        for query_key, dict_set in query_set_dict.items():
            if len(dict_set) < min_cluster_size:
                for entry in dict_set:
                    query_set.remove(entry)
                to_remove.add(query_key)
    
    for query_key in to_remove:
        del query_set_dict[query_key]

    db_set_dict = defaultdict(set)
    print("Reading database...", file=sys.stderr)
    with open(database, "r") as f2:
        for line in f2:
            split_line = line.rstrip().split("\t")
            query_id = split_line[1]
            rep_id = split_line[0]

            if query_id in query_set:
                db_set_dict[rep_id].add(query_id)

    del query_set

    with open(outpref + "_query_dict.pkl", "wb") as f:
        pickle.dump(query_set_dict, f)

    with open(outpref + "_db_dict.pkl", "wb") as f:
        pickle.dump(db_set_dict, f)

    # compare all sets between query and db dictionaries
    identical_sets = []
    partial_sets = []
    no_match_sets = []
    print("Matching sets...", file=sys.stderr)
    for query_key, query_set in tqdm(query_set_dict).items():
        match_found = False
        for db_key, db_set in db_set_dict.items():
            
            # if found identical set, break and move on to next query
            if query_set == db_set:
                identical_sets.append((query_key, db_key, len(query_set), len(db_set), 1.0, 1.0))
                match_found = True
                break
            else:
                # Calculate the intersection of the two sets
                intersection = query_set & db_set
                
                # Check if the intersection is non-empty and not equal to either set
                if intersection and intersection != query_set and intersection != db_set:
                    intersection_prop_query = len(intersection) / len(query_set)
                    intersection_prop_db = len(intersection) / len(db_set)
                    partial_sets.append((query_key, db_key, len(query_set), len(db_set), intersection_prop_query, intersection_prop_db))
                    match_found = True
        
        if match_found == False:
            no_match_sets.append((query_key, "NA", len(query_set), 0, 0.0, 0.0))
                
    print("Writing files...", file=sys.stderr)

    del query_set_dict
    del db_set_dict

    with open(outpref + "_indentical_matches.txt", "w") as o1:
        o1.write("Query\tDatabase\tQuery_set_length\tDB_set_length\tIntersection_query\tIntersection_db\n")
        for query_key, db_key, query_set_length, db_set_length, intersection_prop_query, intersection_prop_db in identical_sets:
            o1.write(query_key + "\t" + db_key + "\t" + str(query_set_length) + "\t" + str(db_set_length) + "\t" + str(intersection_prop_query) + "\t" + str(intersection_prop_db) + "\n")
    
    with open(outpref + "_partial_matches.txt", "w") as o2:
        o2.write("Query\tDatabase\tQuery_set_length\tDB_set_length\tIntersection_query\tIntersection_db\n")
        for query_key, db_key, query_set_length, db_set_length, intersection_prop_query, intersection_prop_db in partial_sets:
            o2.write(query_key + "\t" + db_key + "\t" + str(query_set_length) + "\t" + str(db_set_length) + "\t" + str(intersection_prop_query) + "\t" + str(intersection_prop_db) + "\n")
   
    with open(outpref + "_no_matches.txt", "w") as o3:
        o3.write("Query\tDatabase\tQuery_set_length\tDB_set_length\tIntersection_query\tIntersection_db\n")
        for query_key, db_key, query_set_length, db_set_length, intersection_prop_query, intersection_prop_db in no_match_sets:
            o3.write(query_key + "\t" + db_key + "\t" + str(query_set_length) + "\t" + str(db_set_length) + "\t" + str(intersection_prop_query) + "\t" + str(intersection_prop_db) + "\n")
  

if __name__ == "__main__":
    main()