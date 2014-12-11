import unittest

from probe_generator.test.test_constants import GENOME, ANNOTATION
from probe_generator.amino_acid_probe import AminoAcidProbe
from probe_generator.amino_acid_indel_probe import AminoAcidIndelProbe

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


class TestAminoAcidIndelProbe(unittest.TestCase):
    """Test case docstring.

    """
    def setUp(self):
        self.insertion = select_check_codons(
                AminoAcidIndelProbe.explode("GHI2:P1-G2 ins M /9", ANNOTATION),
                "CCC",
                "GGG")
        self.deletion = select_check_codons(
                AminoAcidIndelProbe.explode("GHI2:del P1-G2 /9", ANNOTATION),
                "CCC",
                "GGG")
        self.single_deletion = select_check_codons(
                AminoAcidIndelProbe.explode("GHI2:del P1-P1 /9", ANNOTATION),
                "CCC",
                "CCC")
        self.indel = select_check_codons(
                AminoAcidIndelProbe.explode("GHI2:del P1-G2 ins M/9", ANNOTATION),
                "CCC",
                "GGG")

    def test_insertion_probe(self):
        self.assertEqual(
                self.insertion.sequence(GENOME),
                "cccATGGGG")

    def test_deletion_probe(self):
        self.assertEqual(
                self.deletion.sequence(GENOME),
                "aaaacccaa")

    def test_single_deletion_probe(self):
        self.assertEqual(
                self.single_deletion.sequence(GENOME),
                "aaaaGGGcc")

    def test_indel_probe(self):
        self.assertEqual(
                self.indel.sequence(GENOME),
                "aaaATGccc")

    def test_insertion_probe_str(self):
        self.assertEqual(
                str(self.insertion),
                "GHI2:P1-G2insM(ATG)/9_BAZ2_3:13")

    def test_deletion_probe_str(self):
        self.assertEqual(
                str(self.deletion),
                "GHI2:delP1-G2/9_BAZ2_3:10")

    def test_indel_probe_str(self):
        self.assertEqual(
                str(self.indel),
                "GHI2:delP1-G2insM(ATG)/9_BAZ2_3:10")


def select_reference_codon(probes, codon):
    for probe in probes:
        if probe.variant.reference == codon:
            the_probe = probe
            break
    else:
        raise Exception(
                "No such codon. Found: {}".format(probe.variant.reference))
    return the_probe


def select_check_codons(probes, left_check, right_check):
    for probe in probes:
        checks = [index.reference for index in probe.variant.check_sequences]
        if checks == [left_check, right_check]:
            the_probe = probe
            break
    else:
        raise Exception(
                "No such checks. Found {}".format(
                    probe.variant.check_sequences))
    return the_probe
