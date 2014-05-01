import unittest

from probe_generator.test.test_constants import GENOME

from probe_generator.coordinate_probe import CoordinateProbe
from probe_generator.probe import InvalidStatement


class TestCoordinateStatementParse(unittest.TestCase):
    """Test cases for parsing coordinate statements.

    """
    def setUp(self):
        self.statement = "1:4-2/2:3+3"
        self.probe = CoordinateProbe.from_statement(self.statement)
        self.specification = {
                  "chromosome1": '1',
                  "start1":       3,
                  "end1":         4,
                  "chromosome2":  '2',
                  "start2":       3,
                  "end2":         5,
                  "rc_side_1":    False,
                  "rc_side_2":    False,
                  "comment":      ""
                  }

    def test_get_simple_probe_sequence(self):
        self.assertEqual(
            "gtaag",
            self.probe.sequence(GENOME))
    
    def test_print_simple_probe(self):
        self.assertEqual(
            "1:4/2:3",
            str(self.probe))

    def test_parse_is_whitespace_insensitive(self):
        whitspace_probe = CoordinateProbe.from_statement(
                 " 1:4\t-\n2    /2 : 3+\t\t\t3")
        self.assertEqual(
                whitspace_probe.sequence(GENOME),
                "gtaag")
        
    def test_parse_probe_with_comments(self):
        self.specification["comment"] = "-- I'm a comment!"
        self.assertEqual(
            CoordinateProbe.from_statement(
                self.statement + " -- I'm a comment!")._spec,
            self.specification)
        
    def test_print_probe_with_comments(self):
        self.probe._spec["comment"] = "-- I'm a comment!"
        self.assertEqual(
            "1:4/2:3-- I'm a comment!",
            str(self.probe))

    def test_parse_raises_InvalidFormat_on_nonsense_statement(self):
        message = "could not parse coordinate statement 'banana'"
        with self.assertRaisesRegex(InvalidStatement, message):
            CoordinateProbe.from_statement('banana')

    def test_side_1_is_reverse_complemented_when_first_side_is_plus(self):
        self.assertTrue(
                CoordinateProbe.from_statement(
                    "1:100+10/2:200+20")._spec['rc_side_1'])

    def test_side_1_is_not_reverse_complemented_when_first_side_is_minus(self):
        self.assertFalse(
                CoordinateProbe.from_statement(
                    "1:100-10/2:200+20")._spec['rc_side_1'])

    def test_side_2_is_reverse_complemented_when_second_side_is_minus(self):
        self.assertTrue(
                CoordinateProbe.from_statement(
                    "1:100-10/2:200-20")._spec['rc_side_2'])

    def test_side_2_is_not_reverse_complemented_when_second_side_is_plus(self):
        self.assertFalse(
                CoordinateProbe.from_statement(
                    "1:100-10/2:200+20")._spec['rc_side_2'])

    def test_statements_with_unmapped_contigs_parsed_correctly(self):
        try:
            CoordinateProbe.from_statement("GL0021.1:1-25 / GL001234.1:2+25")
        except InvalidStatement:
            self.fail()
