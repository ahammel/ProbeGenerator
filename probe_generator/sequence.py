"""Get the sequence range for a probe.

Provides the `sequence_range` function, which returns a dictionary specifying
the genomic location of a probe sequence, given a row of a UCSC gene table and
a fully-realized specification (one without wild-card characters).

"""
import sys

from probe_generator import annotation

WARNING_MESSAGE = (
        "WARNING: probes generated using the '->' syntax may not have the \n"
        "expected value when the end of the first exon is not joined to the \n"
        "start of the second.\n\n"
        "Double-check that your probe statments are specified correctly.\n")


def sequence_range(specification, row_1, row_2):
    """Return the range of base pairs to be extracted from the genome.

    `specification` is a probe specification, such as is returned by
    `probe_statement.parse`. The `specification` must be fully-realized (i.e.,
    no globs except in the 'bases' field). `row_1` and `row_2` are rows from a
    UCSC annotation table.

    Returns a dict in the format:

        {'range_1': (chromosome1, start1, end1),
         'range_2': (chromosome2, start2, end2),
         'reverse': flag}

    Where the chromosome, start, and end specify the genomic locations of the
    probe sequence, and the `flag` is a boolean indicating whether or not the
    second sequence should be reverse-complemented.

    Raises an `InterfaceError` if the `specification` or either of the `rows`
    are improperly formatted.

    Raises a `NoFeatureError` if the `specification` asks for a feature outside
    of the range of the `row`.

    """
    if specification.get('separator') == '/':
        return _positional_sequence_range(
                specification, row_1, row_2)
    elif specification.get('separator') == '->':
        return _read_through_sequence_range(
                specification, row_1, row_2)
    else:
        raise InterfaceError


def _read_through_sequence_range(specification, row_1, row_2):
    """Return the sequence range for a 'read through' probe specification.

    This strategy returns a sequence range for a fusion that is oriented such
    that the two exons form a single transcriptional unit that can potentially
    be transcribed and translated.

    In practice, this is similar to the positional strategy indicated by the
    '/' operator, but with two differences:

        1. The start and end operators are less relevant, as a read-through
           event only makes sense when it covers the end of the first gene and
           the beginning of the second gene.

        2. The two sequences are rearranged such that the end of the first
           gene touches the start of the second.

    For example:
                                                BAR
                                               |=========>
            ..............................................
            ..............................................
            <-------|
                 FOO


            FOO-/BAR+
                    <----====|
            FOO->BAR
                    ====|<----
    """
    _check_read_through_spec(specification)
    if row_1.get('strand') == '+':
        return _positional_sequence_range(
                specification,
                row_1,
                row_2)
    elif row_1.get('strand') == '-':
        return _positional_sequence_range(
                _flip_specification(specification),
                row_2,
                row_1)


def _flip_specification(specification):
    """Return the specification with all of the fields ending in '1' replaced
    with the corresponding field ending in '2' and vice-versa.

    """
    try:
        return dict(
                specification,
                gene1=specification['gene2'],
                gene2=specification['gene1'],
                feature1=specification['feature2'],
                feature2=specification['feature1'],
                side1=specification['side2'],
                side2=specification['side1'],
                bases1=specification['bases2'],
                bases2=specification['bases1'])
    except KeyError as error:
        raise InterfaceError(error)


def _check_read_through_spec(specification):
    """Raises a warning message if the sides of the specification don't make
    sense.

    This is only an issue for probes specified using the read-through syntax.
    In fact, I may make it illegal to specify sides at all for read-through
    statements in a future version.

    """
    try:
        sides_ok = (specification['side1'] == 'end' and
                    specification['side2'] == 'start')
    except KeyError:
        raise InterfaceError
    if not sides_ok:
        print(WARNING_MESSAGE, file=sys.stderr, end="")


def _positional_sequence_range(specification, row_1, row_2):
    """Return the sequence range for a positional probe specification.

    This strategy returns a sequence range based on the sides of the exons
    specified in the probe statement.

    """
    left_chromosome, right_chromosome  = _get_chromosomes(row_1, row_2)
    (left_start,
     left_end,
     right_start,
     right_end) = _get_base_positions(specification, row_1, row_2)

    reverse_complement_flag = _get_rev_comp_flag(
            specification, row_1, row_2)

    return {'chromosome1': left_chromosome,
            'start1':      left_start,
            'end1':        left_end,
            'chromosome2': right_chromosome,
            'start2':      right_start,
            'end2':        right_end,
            'inversion':   reverse_complement_flag}


def _get_chromosomes(*rows):
    """Return the chromosomes of the rows with the 'chr' prefix removed if it
    exists.

    """
    return (row['chrom'].lstrip('chr') for row in rows)


def _get_base_positions(specification, *rows):
    """Yield the genomic coordinates of a probe, given a set of rows.

    Yields the start and end of each row in the order in which they are passed
    to the function.

    """
    for index, row in enumerate(rows, start=1):
        start, end = _get_base_position_per_row(specification, row, index)
        yield start
        yield end


def _get_base_position_per_row(specification, row, index):
    """Return the start and end positions for one side of an event, given a
    row and its index.

    """
    exon_positions = annotation.exons(row)
    try:
        _, which_exon = specification[
                'feature{}'.format(index)]
        bases = specification['bases{}'.format(index)]
        side = specification['side{}'.format(index)]
        strand = row['strand']
    except KeyError as error:
        raise InterfaceError(str(error))

    exon_start, exon_end = _get_exon(exon_positions, which_exon)

    if bases == '*':
        return exon_start, exon_end
    else:
        return _get_base_pair_range(bases, exon_start, exon_end, side, strand)


def _get_base_pair_range(bases, start, end, side, strand):
    """Return the desired sub-range of a genomic feature, given its range, the
    number of base-pairs required, the side from which the sub-range is to be
    extracted, and the strand of the feature.

    In UCSC genome files, the starting base pairs of exons are given from left
    to right across the '+' strand of the chromosome, regardless of the
    orientation of the gene. The locations of the 'start' and the 'end' of an
    exon are switched for a gene on the minus strand.

    """
    if (side == 'start') == (strand == '+'): # <-- see below
        return start, (start + bases - 1)
    else:
        return (end - bases + 1), end
    # The interpretation of the weird-looking conditional pointed out above is
    # this: the _rightmost_ set of base pairs is the _start_ of a feature on
    # the plus strand, or the _end_ of a feature on the minus strand:
    #
    #                 start |                end |
    #                       --------------------->
    #       + .....................................................
    #       - .....................................................
    #                       <---------------------
    #                       ^ end                ^ start


def _get_exon(positions, index):
    """Return the exon at the (1-based) `index` in the `positions` list.

    """
    try:
        return positions[index-1] # zero-indexed list
    except IndexError:
        raise NoFeatureError(
                "specification requires feature 'exon'[{number!s}], "
                "but row specifies only {length} 'exon'(s)".format(
                    number=index,
                    length=len(positions)))


def _get_rev_comp_flag(specification, row_1, row_2):
    """Determine whether the specification represents an inversion event.

    Take a probe specification and two lines from a UCSC genome annotation.

    Returns True if the sequence from the second feature should be
    reverse-complemented in the probe sequence, else False.

    An inversion is indicated in one of two circumstances:

        1. the two features are on the same strand and are fused at the same
           end:


              |---------->                |----------->
            ..................................................
            ..................................................
                         ^                            ^

        2. the two features are on opposite strands and fused at opposite ends:

              |----------->
            ..................................................
            ..................................................
                          ^              <-------------|
                                                       ^

    In either case, the two features must be inverted relative to one another
    in order to bring the 5' and 3' ends of the feature together.

    """
    side1, side2 = specification['side1'], specification['side2']
    strand1, strand2 = row_1['strand'], row_2['strand']
    if side1 == side2 and strand1 == strand2:
        return True
    elif side1 != side2 and strand1 != strand2:
        return True
    else:
        return False


class InterfaceError(Exception):
    """Raised when an object from another module does not provide the expected
    interface.

    """


class NoFeatureError(Exception):
    """Raised when a specification asks for a feature outside of the range of a
    UCSC table row.

    """
