import unittest

from probe_generator.gene_snp_probe import GeneSnpProbe
from probe_generator.test.test_constants import ANNOTATION, GENOME

class TestGeneSnpProbe(unittest.TestCase):
    def setUp(self):
        self.probe, = GeneSnpProbe.explode(
            "ABC: c.2g>t /4",
            ANNOTATION)

        self.other_probe, = GeneSnpProbe.explode(
            "GHI: c.7A>C /50 -- 2nd exon / - strand",
            ANNOTATION)

    def test_gene_snp_probe_spec(self):
        self.assertEqual(
            {"gene":        "ABC",
             "base":        2,
             "reference":   'g',
             "mutation":    't',
             "bases":       4,
             "transcript":  "FOO",
             "chromosome":  "1",
             "index":       3,
             "comment":     ""},
            self.probe._spec)

    def test_gene_snp_probe_string(self):
        self.assertEqual(
            "ABC:c.2g>t/4_FOO_1:3",
            str(self.probe))

    def test_gene_snp_probe_sequence(self):
        self.assertEqual(
            "ctta",
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
             "index":       12,
             "comment":     "-- 2nd exon / - strand"},
            self.other_probe._spec)
