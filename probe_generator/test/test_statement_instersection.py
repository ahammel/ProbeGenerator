"""Test whether it's possible for a string to match more than one statement
regular expression.

"""
import unittest
import re

from greenery import lego

from probe_generator import (coordinate_probe,
                             snp_probe,
                             exon_probe)


class TestRegexIntersection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.snp_fsm = lego.parse(
                deverbosify(
                    snp_probe._SNP_REGEX.pattern))
        cls.exon_fsm = lego.parse(
                deverbosify(
                    exon_probe._PROBE_STATEMENT_REGEX.pattern))
        cls.coord_fsm = lego.parse(
                deverbosify(
                    coordinate_probe._COORDINATE_STATEMENT_REGEX.pattern))

    def test_coord_fsm_does_not_overlap_snp_fsm(self):
        self.assertEqual(
                self.coord_fsm & self.snp_fsm,
                lego.charclass())
        # lego.charclass() == the empty set of strings

    def test_exon_fsm_does_not_overlap_snp_fsm(self):
        self.assertEqual(
                self.exon_fsm & self.snp_fsm,
                lego.charclass())

    def test_exon_fsm_does_not_overlap_coord_fsm(self):
        self.assertEqual(
                self.exon_fsm & self.coord_fsm,
                lego.charclass())


def deverbosify(regex):
    """Return verbose regular expression string with the comments and
    whitespace stripped out.

    Because of the way greenery deals with escapes (i.e., incorrectly), I
    also escape hyphens in character classes even if they're at the end of the
    class, and un-escape the '#' character (the escape character is required
    in Python verbose regexes for a literal '#').

    """
    sub_strings = []
    for line in regex.split("\n"):
        decommented_line = re.sub(r'[^\\]#.*', '', line)
        whitepace_stripped_line = re.sub(r'\s', '', decommented_line)
        escaped_line = re.sub(r'-]', '\-]', whitepace_stripped_line)
        unescaped_line = re.sub(r'\\#', '#', escaped_line)
        sub_strings.append(unescaped_line)
    return ''.join(sub_strings)
