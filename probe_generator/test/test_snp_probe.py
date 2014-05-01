import unittest

from probe_generator.test.test_constants import GENOME
from probe_generator.snp_probe import SnpProbe


class TestSnpProbe(unittest.TestCase):
    def setUp(self):
        self.probe = SnpProbe.from_statement("1:4 t>g /8")

    def test_snp_probe_string(self):
        self.assertEqual(
                "1:4_t>g/8",
                str(self.probe))
        
    def test_snp_probe_string_with_comments(self):
        self.assertEqual(
            "1:4_t>g/8--comment",
            str(SnpProbe.from_statement("1:4 t>g /8 --comment")))

    def test_snp_probe_sequence_simple_test(self):
        self.assertEqual(
                "acggacgt",
                self.probe.sequence(GENOME))

    def test_snp_negative_strand_mutation(self):
        rc_probe = SnpProbe.from_statement("2:4 t>a /8")
        self.assertEqual(
                "aaatgggg",
                rc_probe.sequence(GENOME))

    def test_snp_probe_is_case_insensitive(self):
        upper_probe = SnpProbe.from_statement("1:4 T>G /8")
        self.assertEqual(
                "acgGacgt",
                upper_probe.sequence(GENOME))
