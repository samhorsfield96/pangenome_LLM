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

unique_IDs = set()
sequences = []
with open(input_file_path, 'r') as f:
    for line in f:
        while True:
            line = f.readline()
            if not line:
                break
            gene = line.rstrip()
            gene_ID = line.rstrip().split(" [SEP] ")[0]
            #sequences.append(gene)
            unique_IDs.add(gene_ID)

ID_tokens = [" ".join(unique_IDs) + " A T G C"]
#print(ID_tokens)

tokenizer = Tokenizer(models.WordLevel(unk_token="[UNK]"))
tokenizer.pre_tokenizer = pre_tokenizers.CharDelimiterSplit(" ")
trainer = trainers.WordLevelTrainer(
    special_tokens=["[UNK]", "[CLS]", "[SEP]", "[PAD]", "[MASK]", "[END]"]
)


tokenizer.train_from_iterator(ID_tokens, trainer)
#tokenizer.add_tokens(["A", "C", "G", "T"])
tokenizer.save("tokens.bin")  # Save the trained tokenizer
print("Finished generating tokeniser. Now encoding sequences...")

#print(sequences[0])
#genes = tokenizer.encode_batch(sequences)
#genes = [gene.ids for gene in genes]
#print(genes)

genes = []
with open(input_file_path, 'r') as f:
    for line in f:
        while True:
            line = f.readline()
            if not line:
                break
            gene = line.rstrip()

            encoded = tokenizer.encode(gene).ids
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