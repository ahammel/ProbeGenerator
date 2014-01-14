import unittest
import os

import probe_generator.reference as reference

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
        "/",
        "genesis",
        "scratch",
        "validations",
        "test_probe_generator",
        "test_genome.fa"
        )


class TestReferenceBases(unittest.TestCase):
    """Test cases for the reference.bases function.

    """
    def setUp(self):
        self.ref_genome = MOCK_REFERENCE_GENOME

    def test_bases_returns_inclusive_base_pair_range(self):
        """
        Reference.slice('chr', p, q) should return the base pair sequence
        from the pth base to the qth base inclusive.

        """
        self.assertEqual(
                reference.bases(self.ref_genome, '1', 3, 8),
                'AACCCC')
        self.assertEqual(
                reference.bases(self.ref_genome, 'X', 16, 16),
                't')

    def test_bases_raises_MissingChromosome_when_chromosome_key_missing(self):
        message = "no such chromosome: 'banana'"
        with self.assertRaisesRegex(reference.MissingChromosome, message):
            reference.bases(self.ref_genome, 'banana', 1, 2)

    def test_bases_raises_InvalidRange_on_nonsensical_range(self):
        message = ("unsupported values for `start` and `end`: "
                   "\('foo', 'bar'\)")
        with self.assertRaisesRegex(reference.InvalidRange, message):
            reference.bases(self.ref_genome, '1', 'foo', 'bar')

    def test_bases_raises_NonContainedRange_on_range_outside_of_chromsome(self):
        message = "range \[1:100\] outside the range of chromosome '1'"
        with self.assertRaisesRegex(reference.NonContainedRange, message):
            reference.bases(self.ref_genome, '1', 1, 100)


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
    """Inegration tests for reference.reference_genome and reference.bases.

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

    Calls reference.refernce_genome on production data.

    """
    @classmethod
    def setUpClass(cls):
        with open(PRODUCTION_GENOME_FILE) as handle:
            cls.ref_genome = reference.reference_genome(handle)

    def test_ref_genome_has_proper_chromosomes(self):
        """
        22 autosomes + 2 sex chromosomes

        """
        chromosomes = [
                "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                "12", "13", "14", "15", "16", "17", "18", "19", "20", "21",
                "22", "X", "Y", "MT"]
        for chr in chromosomes:
            self.assertIn(chr, self.ref_genome)

    def test_all_bases_have_only_base_pair_sequences(self):
        for sequence in self.ref_genome.values():
            self.assertNotRegex(
                    sequence,
                    '[^ACGTN]')
