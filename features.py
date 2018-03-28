import numpy as np
import pandas as pd

aa_feats = pd.read_csv("./cpdb2_aa_features.csv", index_col=0)
ss_feats = pd.read_csv("./cpdb2_ss_features.csv", index_col=0)


def prot_to_vector(seq: str, kind: str) -> np.ndarray:
    """Concatenate the amino acid features for each position of the sequence.
    Args:
        seq: A string representing an amino acid sequence.
        kind: One of "aa" or "ss" to indicate amino acid or structure features.

    Returns:
        A numpy array of features, shape (len(seq), features)"""

    # convert to uppercase
    seq = seq.upper()

    if kind == "aa":
        feats = aa_feats
        # replace entries of "B" with "X" in amino acid sequences
        seq = seq.replace("B", "X")
        # replace entries of "J" with "X" in amino acid sequences
        seq = seq.replace("J", "X")
        # replace entries of "O" with "X"
        seq = seq.replace("O", "X")
        # replace entries of "U" with "X"
        seq = seq.replace("U", "X")
        seq = seq.replace("Z", "X")
    elif kind == "ss":
        feats = ss_feats
    else:
        raise ValueError("Invalid argument for argument 'kind'")

    try:
        chain = [feats.loc[pos].values for pos in seq]
    except KeyError as e:
        print(e)
        raise ValueError("Invalid string character encountered in prot_to_vector, kind = %s" % kind)

    return np.concatenate(chain, axis=0).reshape(len(seq), -1)
