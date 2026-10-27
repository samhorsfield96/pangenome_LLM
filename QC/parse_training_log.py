import argparse
import re

def get_options():
    description = "Parses logs from model training"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python parse_training_log.py')
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
    epoch = -1
    with open(infile, "r") as f:
        prev_epoch = -1
        for line in f:
            # ensure on loss line
            if "Loss:" not in line:
                continue

            # format line
            split_line = line.rstrip().split()
            if "Test Loss" in line:
                match = re.search(r'Test Loss:\s*(.*)', line)
                if match:
                    result = match.group(1)
                    split_result = result.rstrip().replace(',', '').split()
                    loss = split_result[0]
                    perplexity = split_result[2]
                    accuracy = split_result[4]
                    precision = split_result[6]
                    recall = split_result[8]
                    F1 = split_result[10]
                    tmp_tuple = (loss, perplexity, accuracy, precision, recall, F1)
                    test_list.append(tmp_tuple)
            else: 
                curr_epoch = split_line[1]
                # check if at new epoch
                if prev_epoch != curr_epoch:
                    epoch += 1
                    prev_epoch = curr_epoch
                if "Training Loss" in line:
                    match = re.search(r'Training Loss:\s*(.*)', line)
                    if match:
                        result = match.group(1)
                        split_result = result.rstrip().replace(',', '').split()
                        #print(result)
                        loss = split_result[0]
                        perplexity = split_result[2]
                        learning_rate = split_result[5]
                        tmp_tuple = (epoch, loss, perplexity, learning_rate)
                        training_list.append(tmp_tuple)
                elif "Validation Loss" in line:
                    match = re.search(r'Validation Loss:\s*(.*)', line)
                    if match:
                        result = match.group(1)
                        split_result = result.rstrip().replace(',', '').split()
                        loss = split_result[0]
                        perplexity = split_result[2]
                        accuracy = split_result[4]
                        precision = split_result[6]
                        recall = split_result[8]
                        F1 = split_result[10]
                        tmp_tuple = (epoch, loss, perplexity, accuracy, precision, recall, F1)
                        validation_list.append(tmp_tuple)

    with open(outfile, "w") as o:
        o.write("Type\tEpoch\tLoss\tPerplexity\tLearning_rate\tAccuracy\tPrecision\tRecall\tF1\n")
        for index, entry in enumerate(validation_list):
            epoch, loss, perplexity, accuracy, precision, recall, F1 = entry
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Validation", epoch, loss, perplexity, training_list[index][3], accuracy, precision, recall, F1))
        for entry in training_list:
            epoch, loss, perplexity, learning_rate = entry
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Training", epoch, loss, perplexity, learning_rate, "NA", "NA", "NA", "NA"))
        for entry in test_list:
            loss, perplexity, accuracy, precision, recall, F1 = entry
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Test", "NA", loss, perplexity, "NA", accuracy, precision, recall, F1))

if __name__ == "__main__":
    main()