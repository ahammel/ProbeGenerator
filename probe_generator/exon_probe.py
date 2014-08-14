"""Implementation of the probe language specification.

"""
import itertools
import re
import sys

from probe_generator import annotation, transcript
from probe_generator.probe import AbstractProbe, InvalidStatement
from probe_generator.sequence_range import SequenceRange

_PROBE_STATEMENT = r"""    # The regex for one side of a probe statement
        \s*                # non-significant whitespace
        ([a-zA-Z0-9_./-]+) # gene name
        \s*\#\s*           # hash-separator
        exon
        \s*\[\s*           # bracket-delimiter
        (\d+|\*)           # exon number
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


class ExonProbe(AbstractProbe):
    """Probe for a mutation event specified as a fusion of two exons.

    """
    _STATEMENT_SKELETON = (
            "{gene1}#exon[{exon1}]{side1}{bases1}"
            "{separator}"
            "{gene2}#exon[{exon2}]{side2}{bases2}"
            "_{breakpoint1}/{breakpoint2}"
            "_{transcript1}_{transcript2}"
            "{comment}")

    def __init__(self, specification, *transcripts):
        self._spec = specification
        self._rows = transcripts

    def __str__(self):
        return self._STATEMENT_SKELETON.format(**self._spec)

    def get_ranges(self):
        """Return the sequence ranges for an exon probe.

        If necessary, each probe half-sequence will be reverse-complemented so
        that the breakpoint is in the centre of the probe. We
        reverse-complement the first half-sequence if it's the start of an exon
        on the plus strand, or the end of an exon on the minus strand. The
        second half-sequence is reverse-complemented if it's the start of an
        exon on the minus strand or the end of an exon on the plus strand.

        If the arrow separator is used, the two sides of the probe will be
        rearranged if necessary so that the reading frames of the two exons
        will be preserved.

        For example:
                                                BAR
                                               |=========>
            ..............................................
            ..............................................
            <-------|
                 FOO


            FOO-/BAR+
                    <----|====
            FOO->BAR
                    ====|<----

        """
        chromosome1 = self._spec["chromosome1"]
        chromosome2 = self._spec["chromosome2"]
        strand1     = self._spec["strand1"]
        strand2     = self._spec["strand2"]
        side1       = self._spec["side1"]
        side2       = self._spec["side2"]
        start1, end1, start2, end2 = self._get_ranges()
        if self._spec['separator'] == '->' and self._spec['strand1'] == '-':
            start1, start2           = start2, start1
            end1, end2               = end2, end1
            chromosome1, chromosome2 = chromosome2, chromosome1
            strand1, strand2         = strand2, strand1
            side1, side2             = side2, side1
        return (
            SequenceRange(
                chromosome1,
                start1,
                end1,
                reverse_complement=(
                    side1 == strand1),
                ),
            SequenceRange(
                chromosome2,
                start2,
                end2,
                reverse_complement=(
                    strand2 != side2),
                ),
            )

    @staticmethod
    def explode(statement, genome_annotation=None):
        """Given an exon probe statement and a genome annotation return all
        possible, unique, unambiguous probes which fit the statement.

        If two or more possible probes have identical coordinates,
        only the first is returned.

        """
        probes = []

        if genome_annotation is None:
            genome_annotation = []
        cached_specifications = set()
        partial_spec = _parse(statement)
        for spec in _expand(partial_spec, genome_annotation):
            spec_hash = _coord_hash(spec)
            if not spec_hash in cached_specifications:
                cached_specifications.add(spec_hash)
                breakpoint1, breakpoint2 = _get_breakpoints(spec)
                probes.append(ExonProbe(dict(spec,
                                             breakpoint1=breakpoint1,
                                             breakpoint2=breakpoint2)))
        return probes

    def _get_ranges(self):
        """Return a four-tuple of the start and end genomic coordinates of the
        two sides of the probe.

        """
        left_range = _get_range(
            self._spec["exon_range_1"],
            self._spec["side1"],
            self._spec["strand1"],
            self._spec["bases1"])
        right_range = _get_range(
            self._spec["exon_range_2"],
            self._spec["side2"],
            self._spec["strand2"],
            self._spec["bases2"])
        return left_range + right_range


def _parse(probe_statement):
    """Return a partial exon probe specification given a statement in probe
    language.

    Parses the `probe_statement` string, returning a dictionary specifying the
    location and breakpoints of the fusion event specified. The returned
    dictionary has eight fields, as follows:

    {
      'gene(1|2)':     The name of the gene
                       'str'

      'feature(1|2)':  The name and number of the feature e.g.: ('exon', 2)
                       ('str' or '*', int or '*')
                       Currently only the value 'exon' is supported.

      'side(1|2)':     The end of the feature from which to construct the probe
                       'start', 'end', or '*'

      'bases(1|2)':    The number of bases to include in the probe
                       int or '*'
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
     exon_number_1,
     side_1,
     bases_1,
     separator,
     gene_2,
     exon_number_2,
     side_2,
     bases_2,
     comment) = match.groups()
    return {
            'gene1':       gene_1,
            'exon1':       _maybe_int(exon_number_1),
            'side1' :      side_1,
            'bases1':      _maybe_int(bases_1),
            'gene2':       gene_2,
            'exon2':       _maybe_int(exon_number_2),
            'side2' :      side_2,
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
                len(left.exons()),
                len(left.exons()))
        for unglobbed_spec in unglobbed_specs:
            try:
                yield _expand_partial_spec(unglobbed_spec, left, right)
            except transcript.NoFeature as error:
                print("Warning: {!s}".format(error), file=sys.stderr)


def _expand_globs(specification, left_exons=None, right_exons=None):
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

    Raises an ExpandError if the number of features is unspecified
    when the feature number is globbed.

    """
    fields = [] # Globbed fields
    values = [] # Possible values for globbed fields

    if specification['side1'] == '*':
        fields.append('side1')
        values.append(('+', '-'))
    if specification['side2'] == '*':
        fields.append('side2')
        values.append(('+', '-'))

    exon_number_1 = specification['exon1']
    exon_number_2 = specification['exon2']

    if exon_number_1 == '*':
        if left_exons is None:
            raise ExpandError(
                    "number of exons must be specified "
                    "when exon number is globbed")
        fields.append('exon1')
        values.append(n+1 for n in range(left_exons))
    if exon_number_2 == '*':
        if right_exons is None:
            raise ExpandError(
                    "number of exons must be specified "
                    "when exon number is globbed")
        fields.append('exon2')
        values.append(n+1 for n in range(right_exons))

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
    first_exon  = specification['exon1']
    second_exon = specification['exon2']
    return dict(specification,
                transcript1=row_1.name,
                transcript2=row_2.name,
                strand1='+' if row_1.plus_strand else '-',
                strand2='+' if row_2.plus_strand else '-',
                exon_range_1=row_1.exon(first_exon),
                exon_range_2=row_2.exon(second_exon),
                chromosome1=row_1.chromosome,
                chromosome2=row_2.chromosome,
                )


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
        spec['chromosome1'], spec['exon_range_1'], spec['side1'],
        spec['chromosome2'], spec['exon_range_2'], spec['side2']]))


def _get_breakpoints(spec):
    """Return the breakpoint strings ("chromosome:index") of a probe given a
    specification.

    """
    chromosome1, start1, end1, _, _ = spec['exon_range_1']
    chromosome2, start2, end2, _, _ = spec['exon_range_2']
    index1 = start1 - 1 if spec['side1'] == spec['strand1'] else end1
    index2 = start2     if spec['side2'] == spec['strand2'] else end2 - 1
    return ("{}:{}".format(chromosome1, index1),
            "{}:{}".format(chromosome2, index2))


def _get_range(exon_range, side, strand, bases):
    """Return the start and end of the base pair range, given the exon, the
    start of the exon requested, and the number of base pairs.

    """
    chromosome, start, end, _, _ = exon_range
    if (side == '+') == (strand == '+'):
        return (start, start + bases)
    else:
        return (end - bases, end)


class ExpandError(Exception):
    """Raised when `expand` has insufficient information to complete.

    """
