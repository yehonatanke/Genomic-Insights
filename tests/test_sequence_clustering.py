import unittest
from src.analysis.sequence_clustering import cluster_sequences, _get_gc_content


class TestSequenceClustering(unittest.TestCase):

    def test_cluster_sequences(self):
        sequences = {
            "seq1": "ATCG",
            "seq2": "GGCC",
            "seq3": "AATT"
        }
        clusters, features = cluster_sequences(sequences, n_clusters=2)

        self.assertEqual(len(clusters), 3)
        self.assertEqual(features.shape, (3, 2))

        # Check if sequences with similar GC content are in the same cluster
        self.assertEqual(clusters["seq1"], clusters["seq3"])
        self.assertNotEqual(clusters["seq1"], clusters["seq2"])

    def test_get_gc_content(self):
        self.assertEqual(_get_gc_content("ATCG"), 50.0)
        self.assertEqual(_get_gc_content("GGCC"), 100.0)
        self.assertEqual(_get_gc_content("AATT"), 0.0)


