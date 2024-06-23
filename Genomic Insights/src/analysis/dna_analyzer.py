import re
from collections import Counter
from Bio import SeqIO
from src.utils.logger import log_message


class DNAAnalyzer:
    def __init__(self, file_path):
        """
        Initialize the DNAAnalyzer with a FASTA file.

        Args:
            file_path (str): Path to the FASTA file.
        """
        self.file_path = file_path
        self.sequences = self._read_fasta()
        log_message("DNAAnalyzer initialized with file: " + file_path)

    def _read_fasta(self):
        """
        Read the FASTA file and return a dictionary of sequences.

        Returns:
            dict: A dictionary with sequence IDs as keys and sequences as values.
        """
        sequences = {record.id: str(record.seq) for record in SeqIO.parse(self.file_path, "fasta")}
        log_message(f"Read {len(sequences)} sequences from FASTA file")
        return sequences

    def get_gc_content(self, sequence):
        """
        Calculate the GC content of a given sequence.

        Args:
            sequence (str): DNA sequence.

        Returns:
            float: GC content percentage.
        """
        gc_count = sequence.count('G') + sequence.count('C')
        gc_content = gc_count / len(sequence) * 100
        log_message(f"Calculated GC content: {gc_content:.2f}%")
        return gc_content

    def find_motifs(self, sequence, motif):
        """
        Find all occurrences of a motif in a sequence.

        Args:
            sequence (str): DNA sequence.
            motif (str): Motif to search for.

        Returns:
            list: List of starting positions of the motif.
        """
        positions = [m.start() for m in re.finditer(f'(?={motif})', sequence)]
        log_message(f"Found {len(positions)} occurrences of motif '{motif}'")
        return positions

    def get_kmer_frequency(self, sequence, k=3):
        """
        Calculate the frequency of k-mers in a sequence.

        Args:
            sequence (str): DNA sequence.
            k (int): Length of k-mer.

        Returns:
            Counter: A Counter object with k-mers and their frequencies.
        """
        kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
        frequency = Counter(kmers)
        log_message(f"Calculated frequency of {k}-mers")
        return frequency
