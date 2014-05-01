import unittest
import re

from probe_generator.exon_probe import ExonProbe
from probe_generator.probe import InvalidStatement


class AbstractExonStatementTestCase(unittest.TestCase):
    """Provides setUp function for testcases in the `probe_statement` module.

    """
    def setUp(self):
        self.probe_statement = (
                "ABC#exon[1]-2 / DEF#exon[3]+3")
        self.probe_specification = {
                    'gene1':       'ABC',
                    'feature1':    ('exon', 1),
                    'side1':       'end',
                    'bases1':      20,
                    'gene2':       'DEF',
                    'feature2':    ('exon', 3),
                    'side2':       'start',
                    'bases2':      30,
                    'separator':   '/',
                    'chromosome1': '1',
                    'breakpoint1': 4,
                    'chromosome2': '2',
                    'breakpoint2': 1,
                }
        self.annotation = [
                {'name':       'FOO',
                 'proteinID':  'ABC',
                 'exonStarts': '1,3,',
                 'exonEnds':   '2,4,',
                 'strand':     '+',
                 'chrom':      '1'},
                {'name':       'BAR',
                 'proteinID':  'DEF',
                 'exonStarts': '7,9,11,',
                 'exonEnds':   '8,10,12,',
                 'strand':     '+',
                 'chrom':      '2'}]



class TestExonStatementParsing(AbstractExonStatementTestCase):
    """Test cases for probe language parsing functionality

    """
    def test_nonsense_probe_statement_raises_exception(self):
        with self.assertRaisesRegex(InvalidStatement, "banana"):
            for _ in ExonProbe.explode("banana"):
                pass

    def test_partial_probe_statement_raises_InvalidStatement(self):
        with self.assertRaises(InvalidStatement):
            for _ in ExonProbe.explode(self.probe_statement[:10]):
                pass

    def test_alphanumerics_and_punctuation_are_valid_in_gene_names(self):
        """
        Certain punctuation characters as well as alphanumerics are valid in
        the names of genes.

        """
        try:
            statement = "abc123#exon[1] -1 / b.a/n_a-na#exon[2] -3"
            for _ in ExonProbe.explode(statement):
                pass
        except InvalidStatement:
            self.fail("Statement could not be parsed")

    def test_parse_probe_statement_with_read_through_separator(self):
        try:
            statement = "ABC#exon[1] -20 -> DEF#exon[3] +30"
            for _ in ExonProbe.explode(statement):
                pass
        except InvalidStatement:
            self.fail("Statement could not be parsed")

    def test_parsing_a_non_exon_feature_raises_not_supported_message(self):
        message = re.escape(
                   "could not parse 'FOO#intron[1]+1/BAR#*[2]-1': "
                   "currently only exons are supported")
        with self.assertRaisesRegex(InvalidStatement, message):
            statement = 'FOO#intron[1]+1/BAR#*[2]-1'
            for _ in ExonProbe.explode(statement):
                pass


class TestExplode(AbstractExonStatementTestCase):
    """Test case for the `ExonProbe` function.

    """
    def test_explode_returns_one_statement_for_fully_realized_statement(self):
        self.assertEqual(
                len(list(
                    ExonProbe.explode(self.probe_statement, self.annotation))),
                1)

    def test_expand_glob_one_side(self):
        statement = "ABC#exon[1]*2 / DEF#exon[3]+3"
        self.assertCountEqual(
                [probe._spec['side1']
                 for probe in ExonProbe.explode(statement, self.annotation)],
                ['start', 'end'])

    def test_expand_glob_both_sides(self):
        statement = "ABC#exon[1]*2 / DEF#exon[3]*3"
        self.assertCountEqual(
                [(probe._spec['side1'], probe._spec['side2'])
                 for probe in ExonProbe.explode(statement, self.annotation)],
                [('start', 'start'), ('start', 'end'),
                 ('end', 'start'),   ('end', 'end')])

    def test_expand_glob_exon_numbers(self):
        statement = "ABC#exon[*]+2 / DEF#exon[*]-3"
        self.assertCountEqual(
                [(probe._spec['feature1'], probe._spec['feature2'])
                  for probe in ExonProbe.explode(statement, self.annotation)],
                 [(('exon', 1), ('exon', 1)),
                  (('exon', 1), ('exon', 2)),
                  (('exon', 1), ('exon', 3)),
                  (('exon', 2), ('exon', 1)),
                  (('exon', 2), ('exon', 2)),
                  (('exon', 2), ('exon', 3))])
