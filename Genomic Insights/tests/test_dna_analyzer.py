import unittest
from src.analysis.dna_analyzer import DNAAnalyzer


class TestDNAAnalyzer(unittest.TestCase):
    def setUp(self):
        self.analyzer = DNAAnalyzer("data/sequences.fasta")

    def test_gc_content(self):
        sequence = "ATGC"
        expected_gc_content = 50.0
        self.assertAlmostEqual(self.analyzer.get_gc_content(sequence), expected_gc_content)

    def test_find_motifs(self):
        sequence = "ATAGCATAGC"
        motif = "ATAGC"
        expected_positions = [0, 5]
        self.assertEqual(self.analyzer.find_motifs(sequence, motif), expected_positions)

    def test_kmer_frequency(self):
        sequence = "ATGATG"
        expected_frequency = {"ATG": 2, "TGA": 1}
        self.assertEqual(dict(self.analyzer.get_kmer_frequency(sequence, k=3)), expected_frequency)


if __name__ == '__main__':
    unittest.main()
