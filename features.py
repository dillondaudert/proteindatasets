import numpy as np
import pandas as pd

def prot_to_vector(aminos_df: pd.DataFrame, seq: str) -> np.ndarray:
    """Concatenate the amino acid features for each position of the sequence.
    Args:
        aminos_df: A pandas DataFrame with the appropriate amino acid features
        seq: A string representing an amino acid sequence.

    Returns:
        A numpy array of features, shape (len(seq), features)"""
    try:
        chain = [aminos_df.loc[residue].values for residue in seq]
    except KeyError:
        raise FeaturesError('Invalid residue encountered in fixed_prot')
    return np.concatenate(chain, axis=0).reshape(len(seq), -1)
