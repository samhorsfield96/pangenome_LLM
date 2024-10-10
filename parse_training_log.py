import argparse

def get_options():
    description = "Compares synteny between simulated and generated genomes"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python synteny_accuracy.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--infile',
                    required=True,
                    help='Infile file name.')
    IO.add_argument('--outfile',
                    default="parsed_log.txt",
                    help='Output file name.')
    return parser.parse_args()

def main():
    options = get_options()
    infile = options.infile
    outfile = options.outfile

    training_list = []
    validation_list = []
    test_list = []
    with open(infile, "r") as f:
        for line in f:
            # format line
            split_line = line.rstrip().replace(',', '').split()
            if "Test" in line:
                loss = split_line[2]
                perplexity = split_line[4]
                accuracy = split_line[6]
                precision = split_line[8]
                recall = split_line[10]
                F1 = split_line[12]
                kappa = split_line[14]
                tmp_tuple = (loss, perplexity, accuracy, precision, recall, F1, kappa)
                test_list.append(tmp_tuple)
            else: 
                epoch = split_line[1]
                loss = split_line[5]
                perplexity = split_line[7]
                if "Training" in line:
                    learning_rate = split_line[10]
                    tmp_tuple = (epoch, loss, perplexity, learning_rate)
                    training_list.append(tmp_tuple)
                elif "Validation" in line:
                    accuracy = split_line[9]
                    precision = split_line[11]
                    recall = split_line[13]
                    F1 = split_line[15]
                    kappa = split_line[17]
                    tmp_tuple = (epoch, loss, perplexity, accuracy, precision, recall, F1, kappa)
                    validation_list.append(tmp_tuple)

    with open(outfile, "w") as o:
        o.write("Type\tEpoch\tLoss\tPerplexity\tLearning_rate\tAccuracy\tPrecision\tRecall\tF1\tKappa\n")
        for index, entry in enumerate(validation_list):
            epoch, loss, perplexity, accuracy, precision, recall, F1, kappa = entry
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Validation", epoch, loss, perplexity, training_list[index][3], accuracy, precision, recall, F1, kappa))
        for entry in training_list:
            epoch, loss, perplexity, learning_rate = entry
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Training", epoch, loss, perplexity, learning_rate, "NA", "NA", "NA", "NA", "NA"))
        for entry in test_list:
            loss, perplexity, accuracy, precision, recall, F1, kappa = entry
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Test", "NA", loss, perplexity, "NA", accuracy, precision, recall, F1, kappa))

if __name__ == "__main__":
    main()