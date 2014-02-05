import unittest

from probe_generator import coordinate_statement as statement


class TestCoordinateStatementParse(unittest.TestCase):
    """Test cases for parsing coordinate statements.

    """
    def setUp(self):
        self.statement = "1:100-10/2:200+20"
        self.specification = {
                  "chromosome1": '1',
                  "start1":       91,
                  "end1":         100,
                  "chromosome2":  '2',
                  "start2":       200,
                  "end2":         219,
                  "inversion":    False}

    def test_parse_simple_coordinate_statement(self):
         self.assertEqual(
                 statement.parse(self.statement),
                 self.specification)

    def test_parse_is_whitespace_insensitive(self):
        self.assertEqual(
                statement.parse(" 1:100\t-\n10    /2 : 200+\t\t\t20"),
                self.specification)

    def test_parse_raises_InvalidFormat_on_nonsense_statement(self):
        message = "could not parse coordinate statement 'banana'"
        with self.assertRaisesRegex(statement.InvalidFormat, message):
            statement.parse('banana')

    def test_parse_returns_an_inversion_when_operations_equal(self):
        self.specification['inversion'] = True
        self.assertEqual(
                statement.parse("1:91+10/2:200+20"),
                self.specification)

    def test_breakpoint_string_returns_string_representation_of_statement(self):
        self.assertEqual(
                statement.breakpoint_string(self.specification),
                "1:100/2:200")
