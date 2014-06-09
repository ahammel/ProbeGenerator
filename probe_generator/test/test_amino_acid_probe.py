import unittest

from probe_generator.test.test_constants import GENOME, ANNOTATION
from probe_generator.amino_acid_probe import AminoAcidProbe

class TestAminoAcidProbe(unittest.TestCase):
    def setUp(self):
        self.probe, = AminoAcidProbe.explode("GHI: P2M /9", ANNOTATION)

    def test_amino_acid_probe_sequence(self):
        self.assertEqual(
            self.probe.sequence(GENOME),
            "ccCATcccc")

    def test_amino_acid_probe_string(self):
        self.assertEqual(str(self.probe), "GHI:P2M(ATG)/9_BAZ_3:15")

    def test_explode_gives_one_probe_per_possible_mutation_codon(self):
        self.assertEqual(
            len(list(AminoAcidProbe.explode("GHI: p2c /9", ANNOTATION))), 2)
        self.assertEqual(
            len(list(AminoAcidProbe.explode("GHI: p2* /9", ANNOTATION))), 3)
        self.assertEqual(
            len(list(AminoAcidProbe.explode("GHI: p2l /9", ANNOTATION))), 6)
