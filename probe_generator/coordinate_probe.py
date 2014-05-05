"""Parse genomic coordinates from simple, human-readable statements.

"""
import re

from probe_generator import reference
from probe_generator.probe import InvalidStatement

_COORDINATE = r"""
    \s*
    ([a-zA-Z0-9.]+) # chromosome
    \s*:\s*         # colon-separator
    (\d+)           # start
    \s*
    ([+-])          # side
    \s*
    (\d+)           # range
    \s*
"""

_COORDINATE_STATEMENT_REGEX = re.compile(r"""
        {0}
        /
        {0}
        (--.*|\s*)       # comment
        """.format(_COORDINATE), re.VERBOSE)

_COORDINATE_STATEMENT_SKELETON = ("{chromosome1}:{end1}/"
                                  "{chromosome2}:{start2}"
                                  "{comment}")


class CoordinateProbe(object):
    """A probe for a fusion event using the coordinates for the breakpoints.

    Coordinate specifications are dictionaries in the format:

        {
        'chromosome(1|2)': str
        'start(1|2)':      int
        'end(1|2)':        int
        'rc_side_(1|2)':   bool
        }

    'start' and 'end' values in the specification are 1-based inclusive ranges.

    """
    def __init__(self, specification):
        self._spec = specification

    def __str__(self):
        return _COORDINATE_STATEMENT_SKELETON.format(**self._spec)

    def sequence(self, genome):
        return reference.bases_from_coordinate(self._spec, genome)

    @staticmethod
    def from_statement(statement):
        """Parse a coordinate probe from a string.

        Statements are in the form:

            "<chr>:<start>(-|+)<range>/<chr>:<start>(-|+)<range>"

        Where <chr> is a chromosome, <start> is the start of the probe, and <range>
        is the number of flanking bases to include in the probe. '+' indicates that
        the flanking bases should be *after* (i.e., the indices of the bases are
        larger than the start) the starting base, and '-' indicates that they
        should be *before*.

        """
        spec = _parse(statement)
        return CoordinateProbe(spec)


def _parse(statement):
    """Return a coordinate specification from a statement.

    """
    match = _COORDINATE_STATEMENT_REGEX.match(statement)
    if not match:
        raise InvalidStatement(
                "could not parse coordinate statement {!r}".format(
                    statement))
    (chr_1,
     start_1,
     operation_1,
     range_1,
     chr_2,
     start_2,
     operation_2,
     range_2,
     comment) = match.groups()

    start_1, end_1 = _parse_range(start_1, operation_1, range_1)
    start_2, end_2 = _parse_range(start_2, operation_2, range_2)
    return {'chromosome1': chr_1,
            'start1':      start_1,
            'end1':        end_1,
            'chromosome2': chr_2,
            'start2':      start_2,
            'end2':        end_2,
            'rc_side_1':   operation_1 == '+', # 'rc' == reverse complement
            'rc_side_2':   operation_2 == '-',
            'comment':     comment}


def _parse_range(start, operation, bases):
    """Given the start, bases, and operation as strings, return the
    bases of the specification.

    """
    if operation == '+':
        return int(start), (int(start) + int(bases) - 1)
    elif operation == '-':
        return (int(start) - int(bases) + 1), int(start)
    else:
        assert False, ("Operation not in '[+-]'.\n"
                       "That's bad. Contact the maintainer.")
