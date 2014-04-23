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

    def tearDown(self):
        sys.stdout.close()
        sys.stdout = self.stdout_backup

    def test_print_fasta_prints_sequences_in_FASTA_format(self):
        print_probes.print_fasta('foo', 'bar')
        self.assertEqual(
                sys.stdout.getvalue(),
                ">foo\nbar\n")
