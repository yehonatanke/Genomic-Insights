import unittest
from unittest.mock import patch
from src.visualization.plots import plot_gc_content, plot_sequence_clusters, plot_kmer_frequency
from collections import Counter
import numpy as np


class TestPlots(unittest.TestCase):

    @patch('matplotlib.pyplot.savefig')
    def test_plot_gc_content(self, mock_savefig):
        sequences = {"seq1": "ATCG", "seq2": "GGCC"}
        gc_contents = [50.0, 100.0]
        plot_gc_content(sequences, gc_contents)
        mock_savefig.assert_called_once_with('gc_content.png')

    @patch('matplotlib.pyplot.savefig')
    def test_plot_sequence_clusters(self, mock_savefig):
        features = np.array([[50.0, 100], [75.0, 200]])
        clusters = np.array([0, 1])
        plot_sequence_clusters(features, clusters)
        mock_savefig.assert_called_once_with('sequence_clusters.png')

    @patch('matplotlib.pyplot.savefig')
    def test_plot_kmer_frequency(self, mock_savefig):
        kmer_freq = Counter({'ATG': 3, 'CGT': 2, 'TTA': 1})
        plot_kmer_frequency(kmer_freq)
        mock_savefig.assert_called_once_with('kmer_frequency.png')

