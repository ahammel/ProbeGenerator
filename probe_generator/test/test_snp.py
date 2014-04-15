import unittest

from probe_generator import snp_statement


class TestSnpProbe(unittest.TestCase):
    def setUp(self):
        self.genome = {
                "1": "acgtacgt",
                "2": "aaaagggg",
                }
        self.probe = snp_statement.SnpProbe("1:4 t>g /8", self.genome)

    def test_snp_probe_string(self):
        self.assertEqual(
                "1:4_t>g/8",
                str(self.probe))

    def test_snp_probe_sequence_simple_test(self):
        self.assertEqual(
                "acggacgt",
                self.probe.sequence())

    def test_snp_probe_reverse_complement_test(self):
        rc_probe = snp_statement.SnpProbe("2:4 t>a /8", self.genome)
        self.assertEqual(
                "aaatgggg",
                rc_probe.sequence())

    def test_snp_probe_is_case_insensitive(self):
        upper_probe = snp_statement.SnpProbe("1:4 T>G /8", self.genome)
        self.assertEqual(
                "acgGacgt",
                upper_probe.sequence())
