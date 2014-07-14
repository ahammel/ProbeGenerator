import unittest
import os

from probe_generator import reference, sequence
from probe_generator.test.test_constants import VALIDATION_DATA_DIR

MOCK_GENOME_FILE = [ # input is any iter of strings
    (">1 CM000663.1 Homo sapiens chromosome 1, "
     "GRCh37 primary reference assembly\n"),
    # 1st line of the Ensembl hg19 reference
    "AAAACCCCGGGGTTTT\n",
    #  4^  8^ 12^ 16^
    (">X CM000685.1 Homo sapiens chromosome X, "
     "GRCh37 primary reference assembly\n"),
    "aaaacccc\n",
    "ggggtttt\n"]

MOCK_REFERENCE_GENOME = {
        '1': "AAAACCCCGGGGTTTT",
        'X': "aaaaccccggggtttt",
        }

# The test production genome file is a file with 1% of the sequence of the
# actual genome reference in production, not including the unmapped contigs
PRODUCTION_GENOME_FILE = os.path.join(
        VALIDATION_DATA_DIR, "test_genome.fa")


class TestReferenceBases(unittest.TestCase):
    """Test cases for the reference.bases function.

    """
    def setUp(self):
        self.ref_genome = MOCK_REFERENCE_GENOME

    def test_bases_returns_base_pair_range(self):
        self.assertEqual(
                reference.bases(
                    sequence.SequenceRange('1', 2, 8),
                    self.ref_genome),
                'AACCCC')
        self.assertEqual(
                reference.bases(
                    sequence.SequenceRange('X', 15, 16),
                    self.ref_genome),
                't')
        self.assertEqual(
                reference.bases(
                    sequence.SequenceRange('1', 0, 16),
                    self.ref_genome),
                'AAAACCCCGGGGTTTT')

    def test_bases_with_reverse_complement(self):
        self.assertEqual(
                reference.bases(
                    sequence.SequenceRange(
                        '1', 2, 8, reverse_complement=True),
                    self.ref_genome),
                'GGGGTT')
        self.assertEqual(
                reference.bases(
                    sequence.SequenceRange(
                        'X', 15, 16, reverse_complement=True),
                    self.ref_genome),
                'a')
        self.assertEqual(
                reference.bases(
                    sequence.SequenceRange(
                        '1', 0, 16, reverse_complement=True),
                    self.ref_genome),
                'AAAACCCCGGGGTTTT')

    def test_bases_raises_MissingChromosome_when_chromosome_key_missing(self):
        message = "no such chromosome: 'banana'"
        with self.assertRaisesRegex(reference.MissingChromosome, message):
            reference.bases(sequence.SequenceRange('banana', 1, 2),
                            self.ref_genome)

    def test_bases_raises_NonContainedRange_on_range_outside_of_chromosome(self):
        message = "range \[1:100\] outside the range of chromosome '1'"
        with self.assertRaisesRegex(reference.NonContainedRange, message):
            reference.bases(sequence.SequenceRange('1', 1, 100),
                            self.ref_genome)


class TestReferenceGenome(unittest.TestCase):
    """Test cases for the reference.reference_genome function.

    """
    def test_reference_returns_ref_genome_dict(self):
        self.assertEqual(
                reference.reference_genome(iter(MOCK_GENOME_FILE)),
                MOCK_REFERENCE_GENOME)

    def test_reference_raises_InvalidGenomeFile_on_empty_input(self):
        message = "genome file empty!"
        with self.assertRaisesRegex(reference.InvalidGenomeFile, message):
            reference.reference_genome(iter([]))

    def test_reference_raises_InvalidGenomeFile_on_nonsense_input(self):
        message = "could not parse input: 'banana'"
        with self.assertRaisesRegex(reference.InvalidGenomeFile, message):
            reference.reference_genome(iter(['banana']))


class TestReferenceGenomeBasesIntegration(TestReferenceBases):
    """Integration tests for reference.reference_genome and reference.bases.

    Calls all the tests of the TestReferenceBases case using a reference
    genome parsed from the MOCK_GENOME_FILE instead of one specified in
    advance.

    """
    def setUp(self):
        self.ref_genome = reference.reference_genome(iter(MOCK_GENOME_FILE))


@unittest.skipIf(not os.path.exists(PRODUCTION_GENOME_FILE),
                 "Production genome file not reachable")
class TestReferenceGenomeValidation(unittest.TestCase):
    """Validation tests for reference.reference_genome.

    Calls `reference.reference_genome` on production data.

    """
    @classmethod
    def setUpClass(cls):
        with open(PRODUCTION_GENOME_FILE) as handle:
            cls.ref_genome = reference.reference_genome(handle)

    def test_ref_genome_has_proper_chromosomes(self):
        chromosomes = [
                "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                "12", "13", "14", "15", "16", "17", "18", "19", "20", "21",
                "22", "X", "Y", "MT"]
        self.assertCountEqual(
                chromosomes,
                self.ref_genome.keys())

    def test_all_bases_have_only_base_pair_sequences(self):
        for sequence in self.ref_genome.values():
            self.assertNotRegex(
                    sequence,
                    '[^ACGTN]')
