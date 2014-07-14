import unittest

from probe_generator.gene_snp_probe import GeneSnpProbe
from probe_generator.sequence import SequenceRange
from probe_generator.test.test_constants import ANNOTATION, GENOME

class TestGeneSnpProbe(unittest.TestCase):
    def setUp(self):
        self.probe, = GeneSnpProbe.explode(
            "ABC: c.1g>t /4",
            ANNOTATION)

        self.other_probe, = GeneSnpProbe.explode(
            "GHI: c.6A>C /50 -- 2nd exon / - strand",
            ANNOTATION)

    def test_gene_snp_probe_spec(self):
        self.assertEqual(
            {"gene":        "ABC",
             "base":        1,
             "reference":   'g',
             "mutation":    't',
             "bases":       4,
             "transcript":  "FOO",
             "chromosome":  "1",
             "index":       SequenceRange('1', 2, 3),
             "index_base":  3,
             "comment":     ""},
            self.probe._spec)

    def test_gene_snp_probe_string(self):
        self.assertEqual(
            "ABC:c.1g>t/4_FOO_1:3",
            str(self.probe))

    def test_gene_snp_probe_sequence(self):
        self.assertEqual(
            "ctta",
            self.probe.sequence(GENOME))

    def test_cross_exon_minus_stand_probe(self):
        self.assertEqual(
            {"gene":        "GHI",
             "base":        6,
             "reference":   'A',
             "mutation":    'G',
             "bases":       50,
             "transcript":  "BAZ",
             "chromosome":  "3",
             "index":       SequenceRange('3', 12, 13),
             "index_base":  13,
             "comment":     "-- 2nd exon / - strand"},
            self.other_probe._spec)
