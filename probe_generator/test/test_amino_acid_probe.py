import unittest

from probe_generator.test.test_constants import GENOME, ANNOTATION
from probe_generator.amino_acid_probe import AminoAcidProbe

class TestAminoAcidProbe(unittest.TestCase):
    def setUp(self):
        self.probe = select_reference_codon(
            AminoAcidProbe.explode("GHI: P2M /9", ANNOTATION),
            "CCC")
        self.transcript_probe = select_reference_codon(
            AminoAcidProbe.explode("MNO: G2M [trans]/9", ANNOTATION),
            "GGG")

    def test_amino_acid_probe_sequence(self):
        self.assertEqual(
            self.probe.sequence(GENOME),
            "cccCATccc")

    def test_transcript_probe_sequence(self):
        self.assertEqual(
            self.transcript_probe.sequence(GENOME),
            "aaaATGaaa")

    def test_amino_acid_probe_sequence_even_number_of_bases(self):
        even_probe = select_reference_codon(
            AminoAcidProbe.explode("GHI: P2M /8", ANNOTATION),
            "CCC")
        self.assertEqual(
            even_probe.sequence(GENOME),
            "ccCATccc")

    def test_amino_acid_probe_string(self):
        self.assertEqual(str(self.probe), "GHI:P2M(CCC>ATG)/9_BAZ_3:13")

    def test_transcript_probe_string(self):
        self.assertEqual(str(self.transcript_probe),
                         "MNO:G2M(GGG>ATG)[trans]/9_FROB_3:13")

    def test_explode_gives_one_probe_per_possible_mutation_codon(self):
        self.assertEqual(
            len(list(AminoAcidProbe.explode("GHI: M2W /9", ANNOTATION))), 1)
        self.assertEqual(
            len(list(AminoAcidProbe.explode("GHI: M2* /9", ANNOTATION))), 3)
        self.assertEqual(
            len(list(AminoAcidProbe.explode("GHI: L2* /9", ANNOTATION))), 18)
        self.assertEqual(
            len(list(AminoAcidProbe.explode("GHI: M2X /9", ANNOTATION))), 63)


def select_reference_codon(probes, codon):
    for probe in probes:
        if probe._spec['reference'] == codon:
            the_probe = probe
            break
    else:
        raise Exception("No such codon")
    return the_probe
