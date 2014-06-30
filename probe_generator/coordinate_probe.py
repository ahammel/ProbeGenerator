"""Parse genomic coordinates from simple, human-readable statements.

"""
import re

from probe_generator.probe import AbstractProbe, InvalidStatement
from probe_generator.sequence import SequenceRange

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

_STATEMENT_REGEX = re.compile(r"""
        {0}
        /
        {0}
        (--.*|\s*)       # comment
        """.format(_COORDINATE), re.VERBOSE)


class CoordinateProbe(AbstractProbe):
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
    _STATEMENT_SKELETON = ("{chromosome1}:{breakpoint1}/"
                           "{chromosome2}:{breakpoint2}"
                           "{comment}")

    def __init__(self, specification, *rows):
        """At init time, the start, end, and breakpoints are caluclated and
        added to the specification.

        The start and end values are passed internally to SequenceRange
        objects, while the breakpoints are reported in the string value of the
        probe.

        """
        breakpoint1, breakpoint2     = _get_breakpoints(specification)

        specification['breakpoint1'] = breakpoint1
        specification['breakpoint2'] = breakpoint2

        super().__init__(specification, *rows)


    @staticmethod
    def explode(statement):
        specification = _parse(statement)
        return [CoordinateProbe(specification)]

    def get_ranges(self):
        start1, end1 = _parse_range(
            self._spec['index1'],
            self._spec['operation1'],
            self._spec['bases1'])
        start2, end2 = _parse_range(
            self._spec['index2'],
            self._spec['operation2'],
            self._spec['bases2'])
        return (
            SequenceRange(
                self._spec['chromosome1'],
                start1,
                end1,
                reverse_complement= self._spec['rc_side_1']),
            SequenceRange(
                self._spec['chromosome2'],
                start2,
                end2,
                reverse_complement= self._spec['rc_side_2']))

def _parse(statement):
    """Return a coordinate specification from a statement.

    """
    match = _STATEMENT_REGEX.match(statement)
    if not match:
        raise InvalidStatement(
                "could not parse coordinate statement {!r}".format(
                    statement))
    (chr_1,
     start_1,
     operation_1,
     bases_1,
     chr_2,
     start_2,
     operation_2,
     bases_2,
     comment) = match.groups()

    return {'chromosome1': chr_1,
            'index1':      int(start_1),
            'bases1':      int(bases_1),
            'operation1':  operation_1,
            'chromosome2': chr_2,
            'index2':      int(start_2),
            'bases2':      int(bases_2),
            'operation2':  operation_2,
            'rc_side_1':   operation_1 == '+', # 'rc' == reverse complement
            'rc_side_2':   operation_2 == '-',
            'comment':     comment}


def _parse_range(index, operation, bases):
    """Given the start, bases, and operation as strings, return the
    bases of the specification.

    """
    if operation == '+':
        return (index-1, index+bases-1)
    else:
        return (index-bases, index)


def _get_breakpoints(specification):
    """Given a coordinate probe specification, return the breakpoints.

    """
    return (_get_breakpoint(specification['index1'],
                            specification['operation1'],
                            specification['bases1'],
                            is_left=True),
            _get_breakpoint(specification['index2'],
                            specification['operation2'],
                            specification['bases2'],
                            is_left=False))

def _get_breakpoint(index, operation, bases, is_left):
    """Return the breakpoint of a half-probe.

    The breakpoint differs from the sequence range in that it is reported using
    an inclusive range convention, rather thant the left-inclusive
    right-exclusive convention used internally.

    """
    if operation == '+' and is_left:
        return index + bases
    elif operation == '-' and is_left:
        return index
    elif operation == '+' and not is_left:
        return index
    elif operation == '-' and not is_left:
        return index - bases
