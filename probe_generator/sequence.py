"""Get the sequence range for a probe.

Provides the `sequence_range` function, which returns a dictionary specifying
the genomic location of a probe sequence, given a row of a UCSC gene table and
a fully-realized specification (one without wild-card characters).

"""
import probe_generator.annotation as annotation


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

    Where the chormosome, start, and end specify the genomic locations of the
    probe sequence, and the `flag` is a boolean indicating whether or not the
    second sequence should be reverse-complemented.

    Raises an `InterfaceError` if the `specification` or either of the `rows`
    are improperly formatted.

    Raises a `NoFeatureError` if the `specification` asks for a feature outside
    of the range of the the `row`.

    """
    left_chromosome = row_1['chrom'].lstrip('chr')
    right_chromosome = row_2['chrom'].lstrip('chr')

    left_start, left_end = _get_base_positions(
            specification, row_1, 1)
    right_start, right_end = _get_base_positions(
            specification, row_2, 2)

    reverse_complement_flag = _get_rev_comp_flag(
            specification, row_1, row_2)

    return {'chromosome1': left_chromosome,
            'start1':      left_start,
            'end1':        left_end,
            'chromosome2': right_chromosome,
            'start2':      right_start,
            'end2':        right_end,
            'inversion':   reverse_complement_flag}


def _get_base_positions(specification, row, row_number):
    """Return the genomic coordinates of a probe.

    Returns a 2-tuple of integers.

    """
    exon_postitions = annotation.exons(row)
    try:
        feature_type, which_exon = specification[
                'feature{}'.format(row_number)]
        bases = specification['bases{}'.format(row_number)]
        side = specification['side{}'.format(row_number)]
        strand = row['strand']
    except KeyError as error:
        raise InterfaceError(str(error))

    try:
        exon_start, exon_end = exon_postitions[which_exon-1] # zero-indexed list
    except IndexError:
        raise NoFeatureError(
                "specification requires feature {type!r}[{number!s}], "
                "but row specifies only {length} {type!r}(s)".format(
                    type=feature_type,
                    number=which_exon,
                    length=len(exon_postitions)))

    if bases == '*':
        return exon_start, exon_end
    elif (side == 'start') == (strand == '+'): # <-- see below
        return exon_start, (exon_start + bases - 1)
    else:
        return (exon_end - bases + 1), exon_end
    # In UCSC genome files, the starting base pairs of exons are given from
    # left to right across the '+' strand of the chromosome, regardless of the
    # orientation of the gene. The locations of the 'start' and the 'end' of an
    # exon are switched for a gene on the minus strand.
    #
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
