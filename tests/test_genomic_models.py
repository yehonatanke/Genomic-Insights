import unittest
import numpy as np
from src.models.genomic_models import GenomicSequenceClassifier, MotifPredictor


class TestGenomicModels(unittest.TestCase):
    def setUp(self):
        """Set up test data."""
        # Test DNA sequences
        self.test_sequences = [
            "ATCGATCGATCG",
            "GCTAGCTAGCTA",
            "ATATATATATAT",
            "GCGCGCGCGCGC"
        ]
        
        self.test_labels = [0, 1, 0, 1]
        
        self.test_motif_positions = [
            [2, 5], 
            [1, 4],  
            [0, 3],
            [2, 5]   
        ]

    def test_genomic_sequence_classifier_initialization(self):
        """Test initialization of GenomicSequenceClassifier."""
        classifier = GenomicSequenceClassifier()
        self.assertIsNotNone(classifier.model)
        self.assertIsNotNone(classifier.scaler)

    def test_genomic_sequence_classifier_feature_extraction(self):
        """Test feature extraction in GenomicSequenceClassifier."""
        classifier = GenomicSequenceClassifier()
        features = classifier._sequence_to_features(self.test_sequences[0])
        
        self.assertIsInstance(features, np.ndarray)
      
        self.assertEqual(features[0], 0.5) 
        self.assertEqual(features[1], 0.5)  
        
        self.assertGreater(len(features), 2)

    def test_genomic_sequence_classifier_training(self):
        """Test training of GenomicSequenceClassifier."""
        classifier = GenomicSequenceClassifier()
        classifier.train(self.test_sequences, self.test_labels)
        
        self.assertTrue(hasattr(classifier.model, 'classes_'))

    def test_genomic_sequence_classifier_prediction(self):
        """Test prediction of GenomicSequenceClassifier."""
        classifier = GenomicSequenceClassifier()
        classifier.train(self.test_sequences, self.test_labels)
        
        # Test prediction
        predictions = classifier.predict(self.test_sequences)
        
        # Check if predictions are in correct format
        self.assertEqual(len(predictions), len(self.test_sequences))
        self.assertTrue(all(isinstance(pred, (int, np.integer)) for pred in predictions))

    def test_motif_predictor_initialization(self):
        """Test initialization of MotifPredictor."""
        predictor = MotifPredictor()
        self.assertIsNotNone(predictor.model)
        self.assertIsNotNone(predictor.scaler)

    def test_motif_predictor_feature_extraction(self):
        """Test feature extraction in MotifPredictor."""
        predictor = MotifPredictor()
        features = predictor._sequence_to_features(self.test_sequences[0], window_size=4)
        
        self.assertIsInstance(features, np.ndarray)
        
        expected_windows = len(self.test_sequences[0]) - 3  
        self.assertEqual(len(features), expected_windows * 2)  

    def test_motif_predictor_training(self):
        """Test training of MotifPredictor."""
        predictor = MotifPredictor()
        predictor.train(self.test_sequences, self.test_motif_positions)
        
        # Check if model is trained
        self.assertTrue(hasattr(predictor.model, 'classes_'))

    def test_motif_predictor_prediction(self):
        """Test prediction of MotifPredictor."""
        predictor = MotifPredictor()
        predictor.train(self.test_sequences, self.test_motif_positions)
        
        # Test prediction
        predictions = predictor.predict(self.test_sequences)
        
        self.assertEqual(len(predictions), len(self.test_sequences))
        self.assertTrue(all(isinstance(pred, (int, np.integer)) for pred in predictions))

    def test_invalid_sequence_handling(self):
        """Test handling of invalid sequences."""
        classifier = GenomicSequenceClassifier()
        invalid_sequences = ["", "123", "ATCGX"]  # Invalid sequences
        
        with self.assertRaises(Exception):
            classifier._sequence_to_features(invalid_sequences[0])
        
        with self.assertRaises(Exception):
            classifier._sequence_to_features(invalid_sequences[1])
        
        with self.assertRaises(Exception):
            classifier._sequence_to_features(invalid_sequences[2])

    def test_empty_input_handling(self):
        """Test handling of empty input lists."""
        classifier = GenomicSequenceClassifier()
        predictor = MotifPredictor()
        
        with self.assertRaises(Exception):
            classifier.train([], [])
        
        with self.assertRaises(Exception):
            predictor.train([], [])


if __name__ == '__main__':
    unittest.main() 
