# Uses code from PanGPT (https://github.com/mol-evol/panGPT) developed by James McInerney
import torch
from torch.utils.data import Dataset
import random

class GenomeDataset(Dataset):
    """
    GenomeDataset: A custom PyTorch Dataset for preprocessing genomic sequences.

    This class facilitates loading and preprocessing genomic sequences for training
    and evaluating deep learning models. It utilizes a provided tokenizer to
    convert text sequences into numerical representations suitable for model input.

    Args:
        texts (list): A list of text strings representing the genomic sequences.
        tokenizer (transformers.PreTrainedTokenizer): A tokenizer object for
            tokenizing the text sequences.
        max_length (int): The maximum allowed length for the processed sequences
            (after tokenization and padding).

    Attributes:
        tokenizer (transformers.PreTrainedTokenizer): The tokenizer used for
            tokenization.
        texts (list): The list of original text sequences.
        max_length (int): The maximum allowed length for the processed sequences.

    Methods:
        __len__() -> int: Returns the number of samples in the dataset.
        __getitem__(idx) -> torch.tensor: Returns a preprocessed genomic sequence
            (tensor) at the specified index.
    """

    def __init__(self, texts, tokenizer):
        self.tokenizer = tokenizer
        self.texts = [self.tokenizer.encode(text).ids for text in texts]

    def __len__(self):
        return len(self.texts)

    def batch_sample(self, genome_idx, max_length):
        #text = self.texts[genome_idx]
        encoded = self.texts[genome_idx]
        pos_idx = torch.randint(len(encoded), (1,))[0]

        data = encoded[pos_idx : pos_idx + max_length]
        obs = encoded[pos_idx + 1 : pos_idx + max_length + 1]

        # add padding token, use -100 as described here: 
        # https://www.reddit.com/r/MachineLearning/comments/157zkks/d_should_i_mask_padding_tokens_when_finetuning_a/?rdt=53406
        padded_data = data + [-100] * (
            max_length - len(data)
        )

        padded_obs = obs + [-100] * (
            max_length - len(obs)
        )

        return torch.tensor(padded_data), torch.tensor(padded_obs)