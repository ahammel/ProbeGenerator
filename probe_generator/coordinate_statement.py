"""Parse genomic coordinates from simple, human-readable statements.

"""
import re

_COORDINATE = r"""
    \s*
    (\w+)       #chromosome
    \s*:\s*
    (\d+)       #start
    \s*
    ([+-])
    \s*
    (\d+)       # range
    \s*
"""

_COORDINATE_STATEMENT_REGEX = re.compile(r"""
        {0}
        \/
        {0}
        """.format(_COORDINATE), re.VERBOSE)


def parse(statement):
    """Return a coordinate specification from a statement.

    Statements are in the form:

        "<chr>:<start>(-|+)<range>/<chr>:<start>(-|+)<range>"

    Where <chr> is a chromosome, <start> is the start of the probe, and <range>
    is the number of flanking bases to include in the probe. '+' indicates that
    the flanking bases should be *after* the starting base, and '-' indicates
    that they should be *before*.

    Coordinate specifications are dictionaries in the format:

        {'chromosome(1|2): str
         'start(1|2):      int
         'end(1|2):        int
        }

    'start' and 'end' values in the specification are 1-based inclusive ranges.

    Example:

        >>> parse("1:100-10/2:200+20")
        {"chromosome1":  "1",
         "start1":       91,
         "end1":         100,
         "chromosome2":  "2",
         "start2":       200,
         "end2":         219}

    """
    match = _COORDINATE_STATEMENT_REGEX.match(statement)
    if not match:
        raise InvalidFormat(
                "could not parse coordinate statement {!r}".format(
                    statement))
    (chr_1,
     start_1,
     operation_1,
     range_1,
     chr_2,
     start_2,
     operation_2,
     range_2) = match.groups()

    start_1, end_1 = _parse_range(start_1, operation_1, range_1)
    start_2, end_2 = _parse_range(start_2, operation_2, range_2)
    return {'chromosome1': chr_1,
            'start1':      start_1,
            'end1':        end_1,
            'chromosome2': chr_2,
            'start2':      start_2,
            'end2':        end_2,
            'inversion':   operation_1 == operation_2}


def _parse_range(start, operation, bases):
    """Return the bases of the specification, given the start and bases as
    strings and the operation.

    """
    if operation == '+':
        return int(start), (int(start) + int(bases) - 1)
    elif operation == '-':
        return (int(start) - int(bases) + 1), int(start)
    else:
        assert False, ("Operation not in '[+-]'.\n"
                       "That's bad. Contact the maintainer.")


class InvalidFormat(Exception):
    """Raised when a coordinate statement cannot be parsed.

    """