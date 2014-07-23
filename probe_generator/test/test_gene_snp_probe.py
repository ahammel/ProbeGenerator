import unittest

from probe_generator.gene_snp_probe import GeneSnpProbe
from probe_generator.test.test_constants import ANNOTATION, GENOME

class TestGeneSnpProbe(unittest.TestCase):
    def setUp(self):
        self.probe, = GeneSnpProbe.explode(
            "ABC: c.1c>t /4",
            ANNOTATION)
        self.transcript_probe, = GeneSnpProbe.explode(
            "DEF: c.2C>A [trans] /3",
            ANNOTATION)
        self.other_probe, = GeneSnpProbe.explode(
            "GHI: c.6C>A /5 -- 2nd exon / - strand",
            ANNOTATION)

    def test_gene_snp_probe_string(self):
        self.assertEqual(
            "ABC:c.1c>t/4_FOO_1:2",
            str(self.probe))

    def test_transcript_probe_string(self):
        self.assertEqual(
            "DEF:c.2C>A[trans]/3_BAR_2:10",
            str(self.transcript_probe))

    def test_gene_snp_probe_string_with_comments(self):
        self.assertEqual(
            "GHI:c.6C>T/5_BAZ_3:13-- 2nd exon / - strand",
            str(self.other_probe))

    def test_gene_snp_probe_sequence(self):
        self.assertEqual(
            "atgt",
            self.probe.sequence(GENOME))

    def test_cross_exon_minus_strand_sequence(self):
        self.assertEqual(
            "ccTGG",
            self.other_probe.sequence(GENOME))

    def test_transcript_sequence(self):
        self.assertEqual(
            "gAt",
            self.transcript_probe.sequence(GENOME))
