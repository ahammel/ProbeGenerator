import unittest

from probe_generator.exon_probe import ExonProbe
from probe_generator.probe import InvalidStatement
from probe_generator.test.test_constants import ANNOTATION, GENOME


class AbstractExonStatementTestCase(unittest.TestCase):
    """Provides setUp function for testcases in the `probe_statement` module.

    """
    def setUp(self):
        self.probe_statement = "ABC#exon[1]-2 / DEF#exon[3]+3"
        self.probe_string = "ABC#exon[1]-2/DEF#exon[3]+3_1:2/2:11_FOO_BAR"


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


class TestExplode(AbstractExonStatementTestCase):
    """Test case for the `ExonProbe` explode function.

    """
    def test_explode_fully_realized_statement(self):
        try:
            probe, = ExonProbe.explode(self.probe_statement, ANNOTATION)
        except ValueError as error:
            self.fail("Wrong number of probes created: {}".format(error))
        else:
            self.assertEqual(
                self.probe_string,
                str(probe))

    def test_explode_fully_realized_statement_with_comment(self):
        probe, = ExonProbe.explode(
            self.probe_statement + " -- bladow, comments!",
            ANNOTATION)
        self.assertEqual(
            self.probe_string + "-- bladow, comments!",
            #                    ^ Note: no space
            str(probe))

    def test_expand_glob_one_side(self):
        statement = "ABC#exon[1]*2 / DEF#exon[3]+3"
        self.assertCountEqual(
                [probe._spec['side1']
                 for probe in ExonProbe.explode(statement, ANNOTATION)],
                ['+', '-'])

    def test_expand_glob_both_sides(self):
        statement = "ABC#exon[1]*2 / DEF#exon[3]*3"
        self.assertCountEqual(
                [(probe._spec['side1'], probe._spec['side2'])
                 for probe in ExonProbe.explode(statement, ANNOTATION)],
                [('+', '+'), ('+', '-'),
                 ('-', '+'), ('-', '-')])

    def test_expand_glob_exon_numbers(self):
        statement = "ABC#exon[*]+2 / DEF#exon[*]-3"
        self.assertCountEqual(
                [(probe._spec['exon1'], probe._spec['exon2'])
                  for probe in ExonProbe.explode(statement, ANNOTATION)],
                 [(1, 1), (1, 2), (1, 3),
                  (2, 1), (2, 2), (2, 3)])


class TestSequence(unittest.TestCase):
    """Test cases for the sequence functionality of exon probes.

    """
    def assert_sequence(self, statement, sequence):
        """Assert that the Probe.sequence method returns the sequence when
        called with the global constant GENOME.

        """
        try:
            probe, = ExonProbe.explode(statement, ANNOTATION)
        except ValueError as error:
            self.fail("More than one probe created: {}".format(error))
        else:
            self.assertEqual(
                sequence,
                probe.sequence(GENOME))

    def test_solidus_plus_plus(self):
        self.assert_sequence("ABC#exon[1]+2 / GHI#exon[2]+2", "cgCC")

    def test_solidus_minus_plus(self):
        self.assert_sequence("ABC#exon[1]-2 / GHI#exon[2]+2", "acCC")

    def test_solidus_plus_minus(self):
        self.assert_sequence("ABC#exon[1]+2 / GHI#exon[2]-2", "cgcc")

    def test_solidus_minus_minus(self):
        self.assert_sequence("ABC#exon[1]-2 / GHI#exon[2]-2", "accc")

    def test_arrow_ABC_first(self):
        self.assert_sequence("ABC#exon[1]-2 -> GHI#exon[2]+2", "acCC")

    def test_arrow_GHI_first(self):
        self.assert_sequence("GHI#exon[2]-2 -> ABC#exon[1]+2", "cgcc")


class TestBreakpoints(unittest.TestCase):
    """Test cases for the breakpoints of exon probes.

    """
    def assert_string(self, statement, string):
        """Assert that the Probe.__str__ method returns the 'string'.

        """
        try:
            probe, = ExonProbe.explode(statement, ANNOTATION)
        except ValueError as error:
            self.fail("More than one probe created: {}".format(error))
        else:
            self.assertEqual(
                str(probe),
                string)

    def test_solidus_plus_plus(self):
        self.assert_string(
            "JKL#exon[1]+2 / GHI#exon[2]+2",
            "JKL#exon[1]+2/GHI#exon[2]+2_4:99/3:14_QUX_BAZ")

    def test_solidus_minus_plus(self):
        self.assert_string(
            "JKL#exon[1]-2 / GHI#exon[2]+2",
            "JKL#exon[1]-2/GHI#exon[2]+2_4:150/3:14_QUX_BAZ")

    def test_solidus_plus_minus(self):
        self.assert_string(
            "JKL#exon[1]+2 / GHI#exon[2]-2",
            "JKL#exon[1]+2/GHI#exon[2]-2_4:99/3:10_QUX_BAZ")

    def test_solidus_minus_minus(self):
        self.assert_string(
            "JKL#exon[1]-2 / GHI#exon[2]-2",
            "JKL#exon[1]-2/GHI#exon[2]-2_4:150/3:10_QUX_BAZ")

    def test_arrow_JKL_first(self):
        self.assert_string(
            "JKL#exon[1]-2 -> GHI#exon[2]+2",
            "JKL#exon[1]-2->GHI#exon[2]+2_4:150/3:14_QUX_BAZ")

    def test_arrow_GHI_first(self):
        self.assert_string(
            "GHI#exon[2]-2 -> JKL#exon[1]+2",
            "GHI#exon[2]-2->JKL#exon[1]+2_3:9/4:100_BAZ_QUX")
