import unittest

from probe_generator.test.test_constants import GENOME
from probe_generator.snp_probe import SnpProbe


class TestSnpProbe(unittest.TestCase):
    def setUp(self):
        self.probe, = SnpProbe.explode("1:4 t>g /8")
        self.globbed_target_probes = SnpProbe.explode("1:4 t>* /8")
        self.globbed_reference_probes = SnpProbe.explode("1:4 *>g /8")
        self.globbed_both_probes = SnpProbe.explode("1:4 *>* /8")

    def test_snp_probe_spec(self):
        self.assertEqual(
            self.probe._spec,
            {"chromosome": "1",
             "index":      4,
             "reference":  "t",
             "mutation":   "g",
             "bases":      8,
             "comment":    ""})

    def test_snp_probe_string(self):
        self.assertEqual(
                "1:4_t>g/8",
                str(self.probe))

    def test_snp_probe_string_with_comments(self):
        self.assertEqual(
            ["1:4_t>g/8--comment"],
            [str(probe) for probe in SnpProbe.explode("1:4 t>g /8 --comment")])

    def test_snp_probe_sequence_simple_test(self):
        self.assertEqual(
                "acggacgt",
                self.probe.sequence(GENOME))

    # def test_snp_negative_strand_mutation(self):
    #     rc_probe, = SnpProbe.explode("2:4 t>a /8")
    #     self.assertEqual(
    #             "aaatgggg",
    #             rc_probe.sequence(GENOME))

    def test_snp_probe_is_case_insensitive(self):
        upper_probe, = SnpProbe.explode("1:4 T>G /8")
        self.assertEqual(
                "acgGacgt",
                upper_probe.sequence(GENOME))

    def test_globbed_target_probe_strings(self):
        self.assertCountEqual(
            ["1:4_t>A/8",
             "1:4_t>C/8",
             "1:4_t>G/8"],
            [str(probe) for probe in self.globbed_target_probes])

    def test_globbed_reference_probe_strings(self):
        self.assertCountEqual(
            ["1:4_A>g/8",
             "1:4_C>g/8",
             "1:4_T>g/8"],
            [str(probe) for probe in self.globbed_reference_probes])

    def test_globbed_reference_and_mutation_probes(self):
        self.assertCountEqual(
            ["1:4_A>C/8",
             "1:4_A>G/8",
             "1:4_A>T/8",
             "1:4_C>A/8",
             "1:4_C>G/8",
             "1:4_C>T/8",
             "1:4_G>A/8",
             "1:4_G>C/8",
             "1:4_G>T/8",
             "1:4_T>A/8",
             "1:4_T>C/8",
             "1:4_T>G/8",
             ],
            [str(probe) for probe in self.globbed_both_probes])
