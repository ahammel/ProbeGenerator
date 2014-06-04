import unittest

from probe_generator.gene_snp_probe import GeneSnpProbe
from probe_generator.test.test_constants import ANNOTATION, GENOME

class TestGeneSnpProbe(unittest.TestCase):
    def setUp(self, ):
        self.probe, = GeneSnpProbe.explode(
            "ABC: c.4t>g /8",
            ANNOTATION)

        self.other_probe, = GeneSnpProbe.explode(
            "GHI: c.7A>C /50 -- 2nd exon / - strand",
            ANNOTATION)

    def test_gene_snp_probe_spec(self):
        self.assertEqual(
            {"gene":        "ABC",
             "base":        4,
             "reference":   't',
             "mutation":    'g',
             "bases":       8,
             "transcript":  "FOO",
             "chromosome":  "1",
             "index":       4,
             "comment":     ""},
            self.probe._spec)

    def test_gene_snp_probe_string(self):
        self.assertEqual(
            "ABC:c.4t>g/8_FOO_1:4",
            str(self.probe))

    def test_gene_snp_probe_sequence(self):
        self.assertEqual(
            "acggacgt", # Same as sequence for SNP probe test.
            self.probe.sequence(GENOME))

    def test_cross_exon_minus_stand_probe(self):
        self.assertEqual(
            {"gene":        "GHI",
             "base":        7,
             "reference":   'A',
             "mutation":    'C',
             "bases":       50,
             "transcript":  "BAZ",
             "chromosome":  "3",
             "index":       15,
             "comment":     "-- 2nd exon / - strand"},
            self.other_probe._spec)
