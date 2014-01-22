import unittest
import sys
import io

from probe_generator import print_probes


class TestPrintProbes(unittest.TestCase):
    """Test cases for the `print_probes` module.

    """
    def setUp(self):
        self.annotations = [
                {'name':       'transcript1',
                 'name2':      'FOO',
                 'exonStarts': '10,',
                 'exonEnds':   '20,',
                 'strand':     '+',
                 'chrom':      'chr1'},
                {'name':       'transcript2',
                 'name2':      'BAR',
                 'exonStarts': '30,50',
                 'exonEnds':   '40,60',
                 'strand':     '+',
                 'chrom':      'chr1'}]
        self.statement = "FOO#exon[1]-5/BAR#exon[*]+5"

        self.stdout_backup = sys.stdout
        sys.stdout = io.StringIO()

    def tearDown(self):
        sys.stdout.close()
        sys.stdout = self.stdout_backup

    def test_reverse_complement_returns_reverse_complemented_sequence(self):
        self.assertEqual(
                print_probes.reverse_complement('agctccgaAATTCCGG'),
                'CCGGAATTtcggagct')

    def test_reverse_complement_ingores_Ns(self):
        self.assertEqual(
                print_probes.reverse_complement('aNNNgg'),
                'ccNNNt')

    def test_print_fasta_prints_sequences_in_FASTA_format(self):
        print_probes.print_fasta('foo', 'bar')
        self.assertEqual(
                sys.stdout.getvalue(),
                ">foo\nbar\n")

    def test_explode_statement(self):
        """
        Note that this functions as an integration test when the functions
        which `explode_statement` calls are not mocked in.

        """
        self.assertCountEqual(
                list(print_probes.explode_statements(
                        [self.statement],
                        self.annotations)),
                [({'chromosome1': '1',
                    'start1':      16,
                    'end1':        20,
                    'chromosome2': '1',
                    'start2':      30,
                    'end2':        34,
                    'inversion': False},
                    'FOO#exon[1]-5/BAR#exon[*]+5 transcript1 transcript2'),
                 ({'chromosome1': '1',
                     'start1':      16,
                     'end1':        20,
                     'chromosome2': '1',
                     'start2':      50,
                     'end2':        54,
                     'inversion': False},
                     'FOO#exon[1]-5/BAR#exon[*]+5 transcript1 transcript2')
                  ])
