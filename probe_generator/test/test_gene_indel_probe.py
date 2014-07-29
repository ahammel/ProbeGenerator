import unittest

from probe_generator.gene_indel_probe import GeneIndelProbe
from probe_generator.test.test_constants import ANNOTATION, GENOME


class TestGeneIndelProbe(unittest.TestCase):
    def setUp(self):
        self.ins_probe, = GeneIndelProbe.explode(
            "ABC: c.1 insT /4",
            ANNOTATION)
        self.del_probe, = GeneIndelProbe.explode(
            "ABC: c.1 delC /4",
            ANNOTATION)
        self.indel_probe, = GeneIndelProbe.explode(
            "ABC: c.1 delC insT /4",
            ANNOTATION)
        self.double_indel_probe, = GeneIndelProbe.explode(
            "ABC: c.1 del CG ins \t TA /4",
            ANNOTATION)

    def test_ins_probe_string(self):
        self.assertEqual(
            "ABC:c.1insT/4_FOO_1:2",
            str(self.ins_probe))

    def test_del_probe_string(self):
        self.assertEqual(
            "ABC:c.1delC/4_FOO_1:2",
            str(self.del_probe))

    def test_indel_probe_string(self):
        self.assertEqual(
            "ABC:c.1delCinsT/4_FOO_1:2",
            str(self.indel_probe))

    def test_double_indel_probe_string(self):
        self.assertEqual(
            "ABC:c.1delCGinsTA/4_FOO_1:2",
            str(self.double_indel_probe))

    def test_ins_probe_sequence(self):
        self.assertEqual(
            "aTcg",
            self.ins_probe.sequence(GENOME))

    def test_del_probe_sequence(self):
        self.assertEqual(
            "agta",
            self.del_probe.sequence(GENOME))

    def test_indel_probe_sequence(self):
        self.assertEqual(
            "aTgt",
            self.indel_probe.sequence(GENOME))

    def test_double_indel_probe_sequence(self):
        self.assertEqual(
            "aTAt",
            self.double_indel_probe.sequence(GENOME))