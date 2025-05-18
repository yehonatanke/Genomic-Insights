from src.analysis.dna_analyzer import DNAAnalyzer
from src.analysis.sequence_clustering import cluster_sequences
from src.visualization.plots import plot_gc_content, plot_sequence_clusters, plot_kmer_frequency
from src.utils.logger import setup_logging, log_message


def main():
    setup_logging()
    log_message("Starting genetic analysis")

    # Initialize DNAAnalyzer
    analyzer = DNAAnalyzer("data/sequences.fasta")

    # Analyze GC content
    gc_contents = [analyzer.get_gc_content(seq) for seq in analyzer.sequences.values()]
    plot_gc_content(analyzer.sequences, gc_contents)

    # Find motifs
    motif = "ATAGC"
    for seq_id, sequence in analyzer.sequences.items():
        motif_positions = analyzer.find_motifs(sequence, motif)
        log_message(f"Motif {motif} found in {seq_id} at positions: {motif_positions}")

    # Analyze k-mer frequency
    for seq_id, sequence in analyzer.sequences.items():
        kmer_freq = analyzer.get_kmer_frequency(sequence)
        log_message(f"Analyzed k-mer frequency for {seq_id}")
        plot_kmer_frequency(kmer_freq)

    # Cluster sequences
    clusters, features = cluster_sequences(analyzer.sequences)
    plot_sequence_clusters(features, list(clusters.values()))

    log_message("Genetic analysis completed")


