import os
import pickle
import requests
import numpy as np
from tokenizers import Tokenizer, models, pre_tokenizers, trainers
from sklearn.model_selection import train_test_split
from random import shuffle
from panGPT import GenomeDataset

# download the tiny shakespeare dataset
input_file_path = "grouped_genes.txt"

# Code from PanGPT (https://github.com/mol-evol/panGPT) developed by James McInerney
# Initialize and train the tokenizer using the 'tokenizers' library

tokenizer = Tokenizer(models.WordLevel(unk_token="[UNK]"))
tokenizer.pre_tokenizer = pre_tokenizers.CharDelimiterSplit(" ")
trainer = trainers.WordLevelTrainer(
    special_tokens=["[UNK]", "[CLS]", "[SEP]", "[PAD]", "[MASK]", "[START]", "[END]"], vocab_size=11
)

tokenizer.train_from_iterator([], trainer)
tokenizer.add_tokens(["A", "C", "G", "T"])
tokenizer.save("tokens.bin")  # Save the trained tokenizer

genes = []
with open(input_file_path, 'r') as f:
    while True:
        line = f.readline()
        if not line:
            break
        gene = " ".join(line.rstrip()).replace(",", "[SEP]")

        # add end of contig and assembly tokens
        sequence = "[START] " + gene + " [END]"
        encoded = tokenizer.encode(sequence).ids
        genes.append(encoded)

# Split the data into training and validation sets (80% training, 20% validation)
train_genomes, val_genomes = train_test_split(genes, test_size=0.2, random_state=42)

train_dataset = GenomeDataset(train_genomes, tokenizer, pre_encoded=True)
val_dataset = GenomeDataset(val_genomes, tokenizer, pre_encoded=True)

# save data
with open("train.bin", "wb") as f:
    pickle.dump(train_dataset, f)

with open("val.bin", "wb") as f:
    pickle.dump(val_dataset, f)    

# save the meta information as well, to help us encode/decode later
meta = {
    'vocab_size': tokenizer.get_vocab_size()
}
with open(os.path.join(os.path.dirname(__file__), 'meta.pkl'), 'wb') as f:
    pickle.dump(meta, f)