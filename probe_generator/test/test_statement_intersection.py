"""Test whether it's possible for a string to match more than one statement
regular expression.

"""
import unittest
import re
import itertools

from greenery import lego

import probe_generator

PROBE_MODULES = (
    probe_generator.coordinate_probe,
    probe_generator.snp_probe,
    probe_generator.exon_probe,
    )


class TestRegexIntersection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fsa_list = [
            lego.parse(deverbosify(module._STATEMENT_REGEX.pattern))
            for module in PROBE_MODULES]

    def assert_non_overlapping(self, fsa1, fsa2):
        self.assertEqual(fsa1 & fsa2,
                         lego.charclass()) # <- empty FSA

    def test_statement_regex_mutual_exclusivity(self):
        for fsa1, fsa2 in itertools.combinations(self.fsa_list, 2):
            self.assert_non_overlapping(fsa1, fsa2)


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
