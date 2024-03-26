import os
import pickle
import requests
import numpy as np
from tokenizers import Tokenizer, models, pre_tokenizers, trainers
from sklearn.model_selection import train_test_split
from random import shuffle
from panGPT import GenomeDataset

# Define a custom pre-tokenization rule using a function
def custom_pre_tokenize(text):
    # Split on whitespace and preserve words containing hyphens as single tokens
    return re.findall(r'\b\w+(?:-\w+)+\b|\b\w+\b|\s+', text)

# download the tiny shakespeare dataset
input_file_path = "tokenised_genomes.txt"

unique_chars = set()
genomes = []
with open(input_file_path, 'r') as f:
    for line in f:
        while True:
                line = f.readline()
                if not line:
                    break
                split_line = " " + line.rstrip().split("\t")[1] + " "

                contigs = split_line.split("_")
                # remove empty contigs
                contigs = [x for x in contigs if x != " "]

                # shuffle contig order
                shuffle(contigs)

                contigs = "_".join(contigs)

                # add end of contig and assembly tokens
                sequence = "[START] _" + contigs + "_ [END]"
                genomes.append(sequence)

# Code from PanGPT (https://github.com/mol-evol/panGPT) developed by James McInerney
# Initialize and train the tokenizer using the 'tokenizers' library
unique_tokens = set(token for genome in genomes for token in genome.split())
vocab_size = len(unique_tokens)

tokenizer = Tokenizer(models.WordLevel(unk_token="[UNK]"))
tokenizer.pre_tokenizer = pre_tokenizers.CharDelimiterSplit(" ")
trainer = trainers.WordLevelTrainer(
    special_tokens=["[UNK]", "[CLS]", "[SEP]", "[PAD]", "[MASK]"], vocab_size=vocab_size
)
tokenizer.train_from_iterator(genomes, trainer)
tokenizer.save("tokens.bin")  # Save the trained tokenizer

# Split the data into training and validation sets (80% training, 20% validation)
train_genomes, val_genomes = train_test_split(genomes, test_size=0.2, random_state=42)

train_dataset = GenomeDataset(train_genomes, tokenizer)
val_dataset = GenomeDataset(val_genomes, tokenizer)

# save data
with open("train.bin", "wb") as f:
    pickle.dump(train_dataset, f)

with open("val.bin", "wb") as f:
    pickle.dump(val_dataset, f)    

# save the meta information as well, to help us encode/decode later
meta = {
    'vocab_size': vocab_size
}
with open(os.path.join(os.path.dirname(__file__), 'meta.pkl'), 'wb') as f:
    pickle.dump(meta, f)
