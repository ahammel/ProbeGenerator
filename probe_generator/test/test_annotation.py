import unittest
import os

from probe_generator import annotation
from probe_generator.transcript import Transcript
from probe_generator.sequence import SequenceRange
from probe_generator.test.test_constants import VALIDATION_DATA_DIR, ANNOTATION

MOCK_ANNOTATION_FILE = [ # input is any iterable of strings
        # UCSC annotation files have a header in this format:
        "#foo\tbar\tbaz\n",
        "ding\tdong\tdong\n",
        ]

MOCK_ROW = {'foo': 'ding', 'bar': 'dong', 'baz': 'dong'}

# Validation test data consists of tables downloaded from UCSC with a mock gene
# 'spiked in'
MOCK_REFSEQ_GENES_FILE = os.path.join(
        VALIDATION_DATA_DIR, "test_refseq_genes.txt")
MOCK_UCSC_GENES_FILE = os.path.join(
        VALIDATION_DATA_DIR, "test_ucsc_genes.txt")


class TestTranscript(unittest.TestCase):
    """Test cases for the annotation.exons function

    """
    def setUp(self):
        self.transcript = Transcript(
            {'strand'     : '+',
             'exonStarts' : '3,10,30,50,',
             'exonEnds'   : '5,20,40,60,',
             'chrom'      : '0',
             'name'       : 'FOO',
             'name2'      : 'BAR',
             'cdsStart'   : '11',
             'cdsEnd'     : '59'})

    def test_exons_returns_exon_sequence_ranges(self):
        self.assertEqual(
            self.transcript.exons(),
            [SequenceRange('0', 3, 5),
             SequenceRange('0', 10, 20),
             SequenceRange('0', 30, 40),
             SequenceRange('0', 50, 60)])

    def test_exon_returns_exon_sequence_range_at_one_based_index(self):
        self.assertEqual(
            self.transcript.exon(2),
            SequenceRange('0', 10, 20))

    def test_exons_returns_reversed_positions_when_strand_minus(self):
        self.transcript.plus_strand = False # Note: never do this outside of a
                                             # test.
        self.assertEqual(
            self.transcript.exons(),
            [SequenceRange('0', 50, 60),
             SequenceRange('0', 30, 40),
             SequenceRange('0', 10, 20),
             SequenceRange('0', 3,  5)])

    def test_coding_exons_returns_coding_sequence_ranges(self):
        self.assertEqual(
            self.transcript.coding_exons(),
            [SequenceRange('0', 11, 20),
             SequenceRange('0', 30, 40),
             SequenceRange('0', 50, 59)])


class TestAnnotationValidation(unittest.TestCase):
    """Validation tests for the `annotation` module.

    Validation consists of attempting to retrieve the exons of a mock gene that
    has been manually 'spiked in' to production data. The gene name will always
    be 'MOCK_GENE', and there will always be exactly one exon from bases 2 to 3.

    """
    def assert_mock_gene_in_file(self, annotation_file):
        """Assert that the mock gene is found in the `annotation_file`.

        `annotation_file` is a handle to one of the UCSC annotation files used
        for testing.

        """
        annotations = annotation.parse_ucsc_file(annotation_file)
        matching_features = annotation.lookup_gene("MOCK_GENE", annotations)
        try:
            mock_row, = tuple(matching_features)
            exons = mock_row.exons()
            self.assertEqual(exons, [SequenceRange('0', 2, 3)])
        except ValueError as error:
            self.fail("Unexpected number of mock genes: {}".format(error))

    @unittest.skipIf(not os.path.exists(MOCK_REFSEQ_GENES_FILE),
                     "{} not reachable".format(MOCK_REFSEQ_GENES_FILE))
    def test_refseq_gene_table(self):
        with open(MOCK_REFSEQ_GENES_FILE) as handle:
            self.assert_mock_gene_in_file(handle)

    @unittest.skipIf(not os.path.exists(MOCK_UCSC_GENES_FILE),
                     "{} not reachable".format(MOCK_UCSC_GENES_FILE))
    def test_ucsc_gene_table(self):
        with open(MOCK_UCSC_GENES_FILE) as handle:
            self.assert_mock_gene_in_file(handle)

    def test_nucleotide_index(self):
        transcript1, transcript2, transcript3, _transcript_4 = ANNOTATION
        self.assertEqual(
            transcript1.nucleotide_index(0),
            SequenceRange('1', 0, 1))
        self.assertEqual(
            transcript2.nucleotide_index(1),
            SequenceRange('2', 8, 9))
        transcript_3_indices = [21, 20, 19, 13, 12, 11, 10, 9]
        for base_pair, index in zip(transcript_3_indices, range(10)):
            self.assertEqual(
                transcript3.nucleotide_index(index),
                SequenceRange('3', base_pair, base_pair+1))

    def test_codon_index(self):
        transcript1, transcript2, transcript3, _transcript_4 = ANNOTATION
        self.assertEqual(
            transcript3.codon_index(2), SequenceRange('3', 12, 15)),
        self.assertEqual(
            transcript3.codon_index(3), SequenceRange('3', 9, 12)),
