import unittest

from probe_generator import coordinate_statement as statement


class TestCoordinateStatementParse(unittest.TestCase):
    """Test cases for parsing coordinate statements.

    """
    def setUp(self):
        self.statement = statement.parse("1:100-10/2:200+20")
        self.specification = {
                  "chromosome1": '1',
                  "start1":       90,
                  "end1":         100,
                  "chromosome2":  '2',
                  "start2":       200,
                  "end2":         220}

    def test_parse_simple_coordinate_statement(self):
         self.assertEqual(
                 self.statement,
                 self.specification)

    def test_parse_is_whitespace_insensitive(self):
        self.assertEqual(
                statement.parse(" 1:100\t-\n10    /2 : 200+\t\t\t20"),
                self.specification)

    def test_parse_raises_InvalidFormat_on_nonsense_statement(self):
        message = "could not parse coordinate statement 'banana'"
        with self.assertRaisesRegex(statement.InvalidFormat, message):
            statement.parse('banana')
