import unittest

from probe_generator.transcript import Transcript
from probe_generator.sequence_range import SequenceRange
from probe_generator.test.test_constants import ANNOTATION

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
             SequenceRange('0',  3,  5)])

    def test_coding_exons_returns_coding_sequence_ranges(self):
        self.assertEqual(
            self.transcript.coding_exons(),
            [SequenceRange('0', 11, 20),
             SequenceRange('0', 30, 40),
             SequenceRange('0', 50, 59)])

    def test_nucleotide_index(self):
        transcript1, transcript2, transcript3, *rest = ANNOTATION
        self.assertEqual(
            transcript1.nucleotide_index(1),
            SequenceRange('1', 1, 2))
        self.assertEqual(
            transcript2.nucleotide_index(2),
            SequenceRange('2', 9, 10))
        transcript_3_indices = [22, 21, 20, 14, 13, 12, 11, 10]
        for base_pair, index in zip(transcript_3_indices, range(1,11)):
            self.assertEqual(
                transcript3.nucleotide_index(index),
                SequenceRange('3', base_pair, base_pair+1))

    def test_codon_index(self):
        transcript1, transcript2, transcript3, *rest = ANNOTATION
        self.assertEqual(
            transcript3.codon_index(1),
            SequenceRange('3', 20, 23, reverse_complement=True)),
        self.assertEqual(
            transcript3.codon_index(2),
            SequenceRange('3', 12, 15, reverse_complement=True))

    def test_transcript_range(self):
        self.assertEqual(
            self.transcript.transcript_range(1, 2),
            [SequenceRange('0', 11, 12)])
        self.assertEqual(
            self.transcript.transcript_range(1, 6),
            [SequenceRange('0', 11, 16)])
        self.assertEqual(
            self.transcript.transcript_range(1, 16),
            [SequenceRange('0', 11, 20),
             SequenceRange('0', 30, 36)])

    def test_base_index_nucleotides(self):
        """Assert `base_index` is the inverse of `nucleotide_index`.

        """
        for transcript in ANNOTATION:
            for i in range(1, len(transcript)+1):
                nuc = transcript.nucleotide_index(i)
                base = transcript.base_index(nuc)
                self.assertEqual(i, base, "Gene name: {}".format(
                        transcript.gene_id))
