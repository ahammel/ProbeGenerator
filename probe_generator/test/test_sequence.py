import unittest

from probe_generator import sequence


class TestSequence(unittest.TestCase):
    """Test cases for sequence manipulation functions

    """
    def test_complement_returns_complemented_sequence(self):
        self.assertEqual(
                sequence.complement('agctccgaAATTCCGG'),
                'tcgaggctTTAAGGCC')

    def test_reverse_complement_returns_reverse_complemented_sequence(self):
        self.assertEqual(
                sequence.reverse_complement('agctccgaAATTCCGG'),
                'CCGGAATTtcggagct')

    def test_reverse_complement_ingores_Ns(self):
        self.assertEqual(
                sequence.reverse_complement('aNNNgg'),
                'ccNNNt')


