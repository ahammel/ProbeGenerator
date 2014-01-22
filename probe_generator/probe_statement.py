"""Implementation of the probe language specification.

"""
import re


_PROBE_STATEMENT = r""" # The regex for one side of a probe statement
        \s*             # non-significant whitespace
        ([a-z-_./]+)    # gene name
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

      'side(1|2)':     The end of the feature from which to construct the probe
                       "start", "end", or "*"

      'bases(1|2)':    The number of bases to include in the probe
                       int or "*"
    }

    See the README for the probe language specification.

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
