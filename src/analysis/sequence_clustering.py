import numpy as np
from sklearn.cluster import KMeans
from src.utils.logger import log_message


def cluster_sequences(sequences, n_clusters=3):
    """
    Cluster sequences based on GC content and sequence length.

    Args:
        sequences (dict): Dictionary of sequences.
        n_clusters (int): Number of clusters.

    Returns:
        dict: Dictionary with sequence IDs as keys and cluster labels as values.
        np.array: Array of features used for clustering.
    """
    features = np.array([[_get_gc_content(seq), len(seq)] for seq in sequences.values()])
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(features)

    result = dict(zip(sequences.keys(), clusters))
    log_message(f"Clustered {len(sequences)} sequences into {n_clusters} clusters")
    return result, features


def _get_gc_content(sequence):
    """Helper function to calculate GC content."""
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence) * 100
