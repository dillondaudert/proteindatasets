import numpy as np
import pandas as pd

aa_feats = pd.read_csv("./aa_feats_final.csv", index_col=0)


def prot_to_vector(seq: str) -> np.ndarray:
    """Concatenate the amino acid features for each position of the sequence.
    Args:
        seq: A string representing an amino acid sequence.
    Returns:
        A numpy array of features, shape (len(seq), features)"""

    # convert to uppercase
    seq = seq.upper()

    try:
        chain = [aa_feats.loc[pos].values for pos in seq]
    except KeyError as e:
        print(e)
        raise ValueError("Invalid string character encountered in prot_to_vector")

    return np.concatenate(chain, axis=0).reshape(len(seq), -1)
