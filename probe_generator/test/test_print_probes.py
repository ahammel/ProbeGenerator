import unittest
import sys
import io

from probe_generator import print_probes


class TestPrintProbes(unittest.TestCase):
    """Test cases for the `print_probes` module.

    """
    def setUp(self):
        self.annotations = [
                {'name':       'transcript1', # unique ID from UCSC table
                 'name2':      'FOO',         # gene symbol in RefSeq tables
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

        self.maxDiff = None

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
        left, right = self.annotations
        self.assertEqual(
                list(print_probes.explode_statements(
                    [self.statement],
                    self.annotations)),
                [({'gene1': 'FOO',
                  'feature1': ('exon', 1),
                  'side1': 'end',
                  'bases1': 5,
                  'gene2': 'BAR',
                  'feature2': ('exon', 1),
                  'side2': 'start',
                  'bases2': 5,
                  'separator': '/'},
                  left, right),
                  ({'gene1': 'FOO',
                    'feature1': ('exon', 1),
                    'side1': 'end',
                    'bases1': 5,
                    'gene2': 'BAR',
                    'feature2': ('exon', 2),
                    'side2': 'start',
                    'bases2': 5,
                    'separator': '/'},
                    left, right),
                  ])

    def test_explode_statement_preserves_separator(self):
        exploded_statements = print_probes.explode_statements(
                [self.statement.replace('/', '->')],
                self.annotations)
        self.assertTrue(
                all(spec['separator'] == '->'
                    for spec, _, _ in exploded_statements))
