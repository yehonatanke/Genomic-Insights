import matplotlib.pyplot as plt
import seaborn as sns
from src.utils.logger import log_message


def plot_gc_content(sequences, gc_contents):
    """
    Create a bar plot of GC content for each sequence.

    Args:
        sequences (dict): Dictionary of sequences.
        gc_contents (list): List of GC content values.
    """
    plt.figure(figsize=(12, 6))
    sns.barplot(x=list(sequences.keys()), y=gc_contents)
    plt.title('GC Content in Sequences')
    plt.xlabel('Sequence ID')
    plt.ylabel('GC Content (%)')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('gc_content.png')
    plt.close()
    log_message("Created GC content plot: gc_content.png")


def plot_sequence_clusters(features, clusters):
    """
    Create a scatter plot of sequence clusters.

    Args:
        features (np.array): Array of features used for clustering.
        clusters (np.array): Array of cluster labels.
    """
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(features[:, 0], features[:, 1], c=clusters, cmap='viridis')
    plt.colorbar(scatter)
    plt.title('Sequence Clusters')
    plt.xlabel('GC Content (%)')
    plt.ylabel('Sequence Length')
    plt.savefig('sequence_clusters.png')
    plt.close()
    log_message("Created sequence clusters plot: sequence_clusters.png")


def plot_kmer_frequency(kmer_freq, top_n=10):
    """
    Create a bar plot of top k-mer frequencies.

    Args:
        kmer_freq (Counter): Counter object with k-mer frequencies.
        top_n (int): Number of top k-mers to plot.
    """
    top_kmers = kmer_freq.most_common(top_n)
    kmers, freqs = zip(*top_kmers)

    plt.figure(figsize=(12, 6))
    sns.barplot(x=list(kmers), y=list(freqs))
    plt.title(f'Top {top_n} Most Common k-mers')
    plt.xlabel('k-mer')
    plt.ylabel('Frequency')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('kmer_frequency.png')
    plt.close()
    log_message(f"Created k-mer frequency plot: kmer_frequency.png")
