"""Implementation of the probe language specification.

"""
import re
import itertools


_PROBE_STATEMENT = r""" # The regex for one side of a probe statement
        \s*             # non-significant whitespace
        ([a-z0-9-_./]+) # gene name
        \s*\#\s*        # hash-separator
        ([a-z]+|\*)     # feature name
        \s*\[\s*        # bracket-delimiter
        (\d+|\*)        # the number of the feature
        \s*\]\s*
        ([-+\*])        # first feature side
        \s*
        (\d+|\*)        # number of bases
        \s*
"""

_PROBE_STATEMENT_REGEX = re.compile(r"""
        {0}
        \/ # Separator for the two features
        {0}
        """.format(_PROBE_STATEMENT),
        re.VERBOSE | re.IGNORECASE)


def parse(probe_statement):
    """Return a probe specification given a statement in probe language.

    Parses the `probe_statement` string, returning a dictionary specifying the
    location and breakpoints of the fusion event specified. The returned
    dictionary has eight fields, as follows:

    {
      'gene(1|2)':     The name of the gene
                       "str"

      'feature(1|2)':  The name and number of the feature e.g.: ('exon', 2)
                       ("str" or "*", int or "*")
                       Currently only the value 'exon' is supported.

      'side(1|2)':     The end of the feature from which to construct the probe
                       "start", "end", or "*"

      'bases(1|2)':    The number of bases to include in the probe
                       int or "*"
    }

    See the README for the probe language specification.

    Example:

        >>> parse("FOO#exon[1]+20/BAR#exon[*]-30")
        {'gene1':    'FOO',
         'feature1': ('exon', 1),
         'side1':    'start',
         'bases1':   20,
         'gene2':    'BAR',
         'feature2': ('exon', '*'),
         'bases2':   30,
         'side2':    'end'}

    """
    match = _PROBE_STATEMENT_REGEX.match(probe_statement)
    if not match:
        raise InvalidStatement(
                "Cannot parse probe statement: {!r}".format(
                    probe_statement))
    (gene_1,
     feature_1,
     feature_number_1,
     side_1,
     bases_1,
     gene_2,
     feature_2,
     feature_number_2,
     side_2,
     bases_2) = match.groups()
    feature_1, feature_2 = feature_1.lower(), feature_2.lower()
    if not feature_1 == feature_2 == 'exon':
        raise InvalidStatement("could not parse {!r}: "
                               "currently only exons are supported".format(
                                   probe_statement))
    return {
            'gene1':    gene_1,
            'feature1': (feature_1, _maybe_int(feature_number_1)),
            'side1' :   side_1.replace('+', 'start').replace('-', 'end'),
            'bases1':   _maybe_int(bases_1),
            'gene2':    gene_2,
            'feature2': (feature_2, _maybe_int(feature_number_2)),
            'side2' :   side_2.replace('+', 'start').replace('-', 'end'),
            'bases2':   _maybe_int(bases_2),
            }


def expand(specification, left_features=None, right_features=None):
    """Yield fully-realized probe statements.

    Given a probe_statement with globs, iterate through all the possible
    interpretations of that statement.

    `right_features` and `left_features` are the the total number of features
    expected in specification['feature1'] and specification['feature2'],
    respectively. These arguments must be specified when expaning a
    specification with the feature number globbed.

    Fully realized probe statements may have a glob value for 'bases'.

    If the statement has no globs, return a generator containing only that
    statement.

    WARNING: globbing feature types is not supported. In fact, exons are the
    only feature that are supported.

    """
    fields = [] # Globbed fields
    values = [] # Possible values for globbed fields

    if specification['side1'] == '*':
        fields.append('side1')
        values.append(('start', 'end'))
    if specification['side2'] == '*':
        fields.append('side2')
        values.append(('start', 'end'))

    feature_1, feature_number_1 = specification['feature1']
    feature_2, feature_number_2 = specification['feature2']

    if feature_number_1 == '*':
        if left_features is None:
            raise ExpandError(
                    "number of features must be specified "
                    "when feature number is globbed")
        fields.append('feature1')
        values.append(tuple((feature_1, n+1) for n in range(left_features)))
    if feature_number_2 == '*':
        if right_features is None:
            raise ExpandError(
                    "number of features must be specified "
                    "when feature number is globbed")
        fields.append('feature2')
        values.append(tuple((feature_2, n+1) for n in range(right_features)))

    if fields:
        for parameters in itertools.product(*values):
            new_parameters = dict(zip(fields, parameters))
            yield dict(specification, **new_parameters)
    else:
        yield specification


def _maybe_int(string):
    """Try to parse the `string` to an `int`. Return the string if this fails.

    """
    try:
        return int(string)
    except ValueError:
        return string


class InvalidStatement(Exception):
    """Raised when a probe statement is poorly formatted.

    """


class ExpandError(Exception):
    """Raised when `expand` has insufficient information to complete.

    """
