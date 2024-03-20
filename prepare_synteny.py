import os
import pickle
import requests
import numpy as np
from tokenizers import Tokenizer, models, pre_tokenizers, trainers
from sklearn.model_selection import train_test_split

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
                split_line = line.rstrip().split("\t")
                # add end of contig and assembly tokens
                sequence = split_line[1]
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

# save data
with open("train.bin", "wb") as f:
    pickle.dump(train_genomes, f)

with open("val.bin", "wb") as f:
    pickle.dump(val_genomes, f)    

""" #data = " ".join(data)
# Create the datasets and data loaders for training and validation
train_dataset = GenomeDataset(train_genomes, tokenizer, 256)
val_dataset = GenomeDataset(val_genomes, tokenizer, 256)

#print(f"length of dataset in characters: {len(data):,}")

# get all the unique characters that occur in this text
chars = sorted(list(unique_chars))
vocab_size = len(chars)
#print("all the unique characters:", ''.join(chars))
print(f"vocab size: {vocab_size:,}")

# create a mapping from characters to integers
stoi = { ch:i for i,ch in enumerate(chars) }
itos = { i:ch for i,ch in enumerate(chars) }
def encode(s):
    return [stoi[c] for c in s] # encoder: take a string, output a list of integers
def decode(l):
    return ''.join([itos[i] for i in l]) # decoder: take a list of integers, output a string

# create the train and test splits
n = len(data)
train_data = data[:int(n*0.9)]
val_data = data[int(n*0.9):]

# encode both to integers
train_ids = encode(train_data)
val_ids = encode(val_data)
print(f"train has {len(train_ids):,} tokens")
print(f"val has {len(val_ids):,} tokens")

# export to bin files
train_ids = np.array(train_ids, dtype=np.uint16)
val_ids = np.array(val_ids, dtype=np.uint16)
train_ids.tofile(os.path.join(os.path.dirname(__file__), 'train.bin'))
val_ids.tofile(os.path.join(os.path.dirname(__file__), 'val.bin')) """

# save the meta information as well, to help us encode/decode later
meta = {
    'vocab_size': vocab_size
}
with open(os.path.join(os.path.dirname(__file__), 'meta.pkl'), 'wb') as f:
    pickle.dump(meta, f)
