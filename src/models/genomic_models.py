import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from Bio import SeqIO
from src.utils.logger import log_message


class GenomicSequenceClassifier:
    def __init__(self):
        """Initialize the genomic sequence classifier."""
        self.model = RandomForestClassifier(n_estimators=100, random_state=42)
        self.scaler = StandardScaler()
        log_message("Initialized GenomicSequenceClassifier")

    def _sequence_to_features(self, sequence):
        """
        Convert DNA sequence to numerical features.
        
        Args:
            sequence (str): DNA sequence
            
        Returns:
            numpy.ndarray: Feature vector
        """
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
        at_content = (sequence.count('A') + sequence.count('T')) / len(sequence)
        
        # K-mer frequencies (k=3)
        kmers = [sequence[i:i+3] for i in range(len(sequence)-2)]
        kmer_freq = {}
        for kmer in kmers:
            kmer_freq[kmer] = kmer_freq.get(kmer, 0) + 1
        
        # Normalize k-mer frequencies
        total_kmers = len(kmers)
        kmer_features = [kmer_freq.get(kmer, 0) / total_kmers for kmer in sorted(set(kmers))]
        
        return np.array([gc_content, at_content] + kmer_features)

    def train(self, sequences, labels):
        """
        Train the classifier on genomic sequences.
        
        Args:
            sequences (list): List of DNA sequences
            labels (list): List of corresponding labels
        """
        # Convert sequences to features
        X = np.array([self._sequence_to_features(seq) for seq in sequences])
        
        X_train, X_test, y_train, y_test = train_test_split(
            X, labels, test_size=0.2, random_state=42
        )
        
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        # Train model
        self.model.fit(X_train_scaled, y_train)
        
        # Calculate and log accuracy
        train_score = self.model.score(X_train_scaled, y_train)
        test_score = self.model.score(X_test_scaled, y_test)
        log_message(f"Training accuracy: {train_score:.3f}")
        log_message(f"Testing accuracy: {test_score:.3f}")

    def predict(self, sequences):
        """
        Predict labels for new sequences.
        
        Args:
            sequences (list): List of DNA sequences to predict
            
        Returns:
            list: Predicted labels
        """
        X = np.array([self._sequence_to_features(seq) for seq in sequences])
        X_scaled = self.scaler.transform(X)
        return self.model.predict(X_scaled)


class MotifPredictor:
    def __init__(self):
        """Initialize the motif predictor."""
        self.model = RandomForestClassifier(n_estimators=100, random_state=42)
        self.scaler = StandardScaler()
        log_message("Initialized MotifPredictor")

    def _sequence_to_features(self, sequence, window_size=50):
        """
        Convert sequence to features for motif prediction.
        
        Args:
            sequence (str): DNA sequence
            window_size (int): Size of sliding window
            
        Returns:
            numpy.ndarray: Feature vector
        """
        features = []
        for i in range(0, len(sequence) - window_size + 1):
            window = sequence[i:i+window_size]
            # Calculate window features
            gc_content = (window.count('G') + window.count('C')) / window_size
            at_content = (window.count('A') + window.count('T')) / window_size
            features.extend([gc_content, at_content])
        
        return np.array(features)

    def train(self, sequences, motif_positions):
        """
        Train the motif predictor.
        
        Args:
            sequences (list): List of DNA sequences
            motif_positions (list): List of motif positions for each sequence
        """
        X = []
        y = []
        
        for seq, positions in zip(sequences, motif_positions):
            features = self._sequence_to_features(seq)
            X.append(features)
            y.extend(positions)
        
        X = np.array(X)
        y = np.array(y)
        
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )
        
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        self.model.fit(X_train_scaled, y_train)
        
        train_score = self.model.score(X_train_scaled, y_train)
        test_score = self.model.score(X_test_scaled, y_test)
        log_message(f"Training accuracy: {train_score:.3f}")
        log_message(f"Testing accuracy: {test_score:.3f}")

    def predict(self, sequences):
        """
        Predict motif positions in sequences.
        
        Args:
            sequences (list): List of DNA sequences to predict
            
        Returns:
            list: Predicted motif positions for each sequence
        """
        X = np.array([self._sequence_to_features(seq) for seq in sequences])
        X_scaled = self.scaler.transform(X)
        return self.model.predict(X_scaled) 
