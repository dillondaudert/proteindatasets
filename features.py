import numpy as np
import pandas as pd

aminos_df = pd.read_csv("./aminos.csv", index_col=0)

def prot_to_vector(seq: str) -> np.ndarray:
    """Concatenate the amino acid features for each position of the sequence.
    Args:
        seq: A string representing an amino acid sequence.

    Returns:
        A numpy array of features, shape (len(seq), features)"""
    try:
        chain = [aminos_df.loc[residue].values for residue in seq]
    except KeyError:
        raise FeaturesError('Invalid residue encountered in fixed_prot')
    return np.concatenate(chain, axis=0).reshape(len(seq), -1)
