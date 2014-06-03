"""Implementation of the probe language specification.

"""
import itertools
import re
import sys

from probe_generator import reference, annotation, exon_coordinate
from probe_generator.probe import InvalidStatement

_PROBE_STATEMENT = r"""    # The regex for one side of a probe statement
        \s*                # non-significant whitespace
        ([a-zA-Z0-9_./-]+) # gene name
        \s*\#\s*           # hash-separator
        ([a-zA-Z]+|\*)     # feature name
        \s*\[\s*           # bracket-delimiter
        (\d+|\*)           # the number of the feature
        \s*\]\s*
        ([*+-])            # feature side
        \s*
        (\d+|\*)           # number of bases
        \s*
"""

_SEPARATOR = r"""
    (/      # standard separator
    |
    ->)     # read-through separator
"""

_STATEMENT_REGEX = re.compile(r"""
        {statement}
        {separator}
        {statement}
        (--.*|\s*)   # comment
        """.format(statement=_PROBE_STATEMENT,
                   separator=_SEPARATOR),
        re.VERBOSE)

_STATEMENT_SKELETON = (
        "{gene1}#{feature1[0]}[{feature1[1]}]{side1}{bases1}"
        "{separator}"
        "{gene2}#{feature2[0]}[{feature2[1]}]{side2}{bases2}"
        "_{chromosome1}:{end1}/{chromosome2}:{start2}"
        "_{transcript1}_{transcript2}"
        "{comment}")


class ExonProbe(object):
    """Probe for a mutation event specified as a fusion of two exons.

    """
    def __init__(self, specification):
        self._spec = specification

    def __str__(self):
        side1 = '+' if self._spec['side1'] == 'start' else '-'
        side2 = '+' if self._spec['side2'] == 'start' else '-'
        return _STATEMENT_SKELETON.format(
                **dict(self._spec,
                       side1=side1,
                       side2=side2))

    def sequence(self, genome):
        """Return the sequence of the probe, given a genome reference.

        """
        return reference.bases_from_coordinate(self._spec, genome)

    @staticmethod
    def explode(statement, genome_annotation=None):
        """Given an exon probe statement and a genome annotation return all
        possible, unique, unambiguous probes which fit the statement.

        If two or more possible probes have identical coordinates,
        only the first is returned.

        """
        if genome_annotation is None:
            genome_annotation = []
        cached_specifications = set()
        partial_spec = _parse(statement)
        for spec in _expand(partial_spec, genome_annotation):
            spec_hash = _coord_hash(spec)
            if not spec_hash in cached_specifications:
                cached_specifications.add(spec_hash)
                yield ExonProbe(spec)


def _parse(probe_statement):
    """Return a partial exon probe specification given a statement in probe
    language.

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

    Raises InvalidStatement if the probe_statement cannot be parsed,
    or if a feature aside from an exon is requested.

    """
    match = _STATEMENT_REGEX.match(probe_statement)
    if not match:
        raise InvalidStatement(
                "Cannot parse probe statement: {!r}".format(
                    probe_statement))
    (gene_1,
     feature_1,
     feature_number_1,
     side_1,
     bases_1,
     separator,
     gene_2,
     feature_2,
     feature_number_2,
     side_2,
     bases_2,
     comment) = match.groups()
    feature_1, feature_2 = feature_1.lower(), feature_2.lower()
    if not feature_1 == feature_2 == 'exon':
        raise InvalidStatement("could not parse {!r}: "
                               "currently only exons are supported".format(
                                   probe_statement))
    return {
            'gene1':       gene_1,
            'feature1':    (feature_1, _maybe_int(feature_number_1)),
            'side1' :      side_1.replace('+', 'start').replace('-', 'end'),
            'bases1':      _maybe_int(bases_1),
            'gene2':       gene_2,
            'feature2':    (feature_2, _maybe_int(feature_number_2)),
            'side2' :      side_2.replace('+', 'start').replace('-', 'end'),
            'bases2':      _maybe_int(bases_2),
            'separator':   separator,
            'comment':     comment,
            }


def _expand(specification, genome_annotation):
    """Populate the specification using the values looked up in the
    genome annotation.

    If expanding the specification asks for a feature which is not in
    the annotation, a warning message is printed to standard error.

    If expanding the specification asks for a feature which is not in
    the annotation, a warning message is printed to standard error.

    """
    left_rows = annotation.lookup_gene(
            specification['gene1'], genome_annotation)
    right_rows = annotation.lookup_gene(
            specification['gene2'], genome_annotation)
    for left, right in itertools.product(left_rows, right_rows):
        unglobbed_specs = _expand_globs(
                specification,
                len(annotation.exons(left)),
                len(annotation.exons(right)))
        for unglobbed_spec in unglobbed_specs:
            try:
                yield _expand_partial_spec(unglobbed_spec, left, right)
            except exon_coordinate.NoFeatureError as error:
                print("Warning: {!s}".format(error), file=sys.stderr)


def _expand_globs(specification, left_features=None, right_features=None):
    """Yield fully-realized probe statements.

    Given a probe_statement with globs, iterate through all the possible
    interpretations of that statement.

    `left_features` and `right_features` are the the total number of features
    expected in specification['feature1'] and specification['feature2'],
    respectively. These arguments must be specified when expanding a
    specification with the feature number globbed.

    Fully realized probe statements may have a glob value for 'bases'.

    If the statement has no globs, return a generator containing only that
    statement.

    WARNING: globbing feature types is not supported. In fact, exons are the
    only feature that are supported.

    Raises an ExpandError if the number of features is unspecified
    when the feature number is globbed.

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


def _expand_partial_spec(specification, row_1, row_2):
    """Given a partial specification with no globs, add the unique transcript
    identifiers and the coordinates of the probe to the specification.

    """
    coordinate = exon_coordinate.sequence_range(specification, row_1, row_2)
    return dict(specification,
                transcript1=row_1['name'],
                transcript2=row_2['name'],
                **coordinate)


def _maybe_int(string):
    """Try to parse the `string` to an `int`. Return the string if this fails.

    """
    try:
        return int(string)
    except ValueError:
        return string


def _coord_hash(spec):
    """Return a unique specification identifier.

    """
    return hash(tuple([
        spec['chromosome1'], spec['start1'], spec['end1'],
        spec['chromosome2'], spec['start2'], spec['end2']]))


class ExpandError(Exception):
    """Raised when `expand` has insufficient information to complete.

    """
