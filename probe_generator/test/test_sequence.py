import unittest

from probe_generator import sequence
from probe_generator.sequence_range import SequenceRange


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

    def test_translate(self):
        self.assertEqual(sequence.translate("AAAACGGCT"), "KTA")

    def test_reverse_translate(self):
        self.assertCountEqual(
                sequence.reverse_translate("GM"),
                ["GGAATG", "GGCATG", "GGGATG", "GGTATG"])


class TestSequenceRange(unittest.TestCase):
    def setUp(self):
        self.range_12 = SequenceRange('0', 1, 2)
        self.range_24 = SequenceRange('0', 2, 4)
        self.range_56 = SequenceRange('0', 5, 6)

    def test_concat(self):
        self.assertEqual(
            self.range_12.concat(self.range_24),
            SequenceRange('0', 1, 4))
        self.assertEqual(
            self.range_24.concat(self.range_12),
            SequenceRange('0', 1, 4))
        self.assertRaises(
            ValueError,
            lambda: self.range_12.concat(self.range_56))

    def test_adjacent(self):
        self.assertTrue(
            self.range_12.adjacent(self.range_24))
        self.assertTrue(
            self.range_24.adjacent(self.range_12))
        self.assertFalse(
            self.range_12.adjacent(self.range_56))

    def test_condense(self):
        self.assertEqual(
            SequenceRange.condense(
                self.range_12,
                self.range_24),
            [SequenceRange('0', 1, 4)])
        self.assertEqual(
            SequenceRange.condense(
                self.range_12,
                self.range_24,
                self.range_56),
            [SequenceRange('0', 1, 4),
             self.range_56])

    def test_between(self):
        self.assertEqual(
                SequenceRange.between(
                    self.range_12,
                    self.range_24),
                SequenceRange('0', 2, 2))

    def test_span(self):
        self.assertEqual(
                SequenceRange.span(
                    self.range_12,
                    self.range_24),
                SequenceRange('0', 1, 4))
        self.assertEqual(
                SequenceRange.span(
                    self.range_12,
                    self.range_56),
                SequenceRange('0', 1, 6))

