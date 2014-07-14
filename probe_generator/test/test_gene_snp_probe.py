import unittest

from probe_generator.gene_snp_probe import GeneSnpProbe
from probe_generator.test.test_constants import ANNOTATION, GENOME

class TestGeneSnpProbe(unittest.TestCase):
    def setUp(self):
        self.probe, = GeneSnpProbe.explode(
            "ABC: c.1c>t /4",
            ANNOTATION)

        self.other_probe, = GeneSnpProbe.explode(
            "GHI: c.6C>A /5 -- 2nd exon / - strand",
            ANNOTATION)

    def test_gene_snp_probe_string(self):
        self.assertEqual(
            "ABC:c.1c>t/4_FOO_1:2",
            str(self.probe))

    def test_gene_snp_probe_sequence(self):
        self.assertEqual(
            "atgt",
            self.probe.sequence(GENOME))

    def test_cross_exon_minus_strand_sequence(self):
        self.assertEqual(
            "ccTGG",
            self.other_probe.sequence(GENOME))
