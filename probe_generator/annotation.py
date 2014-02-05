"""Read UCSC annotation tables.

Annotation files can be downloaded from the UCSC table browser:

    http://genome.ucsc.edu/cgi-bin/hgTables

Make sure to use the output format 'all fields from selected table'.

Currently supported tables are:
    RefSeq Genes
    UCSC Genes

Tables can (currently) only be supported if they have the 'exonStarts' and
'exonEnds' fields, and a field identifying the gene name.

To add support for a new table:
    1. Download a copy of the table from UCSC and place it in the validation
       data directory (this directory is listed in the test_constants module).
    2. 'Spike in' a test gene to the data. This gene must have the
       identifier 'MOCK_GENE', the 'exonStarts' value '2,' and the 'exonEnds'
       value '3,'. The values of the other fields are not important. Make sure
       the new row follows the format of the rest of the gene table
       (tab-delimited, correct number of fields).
    3. Add a validation test to the test_annotation.TestAnnotationValidation
       class, following the pattern of the other tests in that class. Make sure
       this test fails properly.
    4. Add the gene identifier field to the _GENE_NAME_FIELDS module constant
       below. Make sure the validation test passes.
    5. Add the name of the gene table to list of supported tables above.

"""
import csv

_GENE_NAME_FIELDS = (
        # the names of the fields which might contain the name of a gene in any
        # of the UCSC file formats listed in the docstring
        'name2',
        'proteinID',
        )


def parse_ucsc_file(handle):
    """Return a csv.DictReader object relating fields of UCSC annotation files
    to values.

    `handle` is an iterator of the lines of file, which must be in the standard
    UCSC format (i.e., tab-delimited file, first line starts with '#' and
    specifies field names.

    """
    lines = (line.lstrip('#') for line in handle)
    return csv.DictReader(lines, dialect='excel-tab')


def lookup_gene(gene_name, ucsc_file):
    """Yield data for features in a `uscs_file` for a specific gene.

    `ucsc_file` is an iterator of dictionaries giving the data from a UCSC gene
    file, as might be returned by `parse_ucsc_file`. Currently supported
    formats are given in the docstring of the `annotation` module.

    """
    for row in ucsc_file:
        if any(row.get(field) == gene_name for field in _GENE_NAME_FIELDS):
            yield row


def exons(row):
    """Return the exon positions of a UCSC annotation feature.

    In a UCSC annotation file, the positions of the starts and ends of exons
    are stored as comma-separated strings:

        '20,30,40,'

    Given a dictionary with this data, we return of a list of tuples:

        (exonStart, exonEnd)

    If the 'strand' of the row is '-', the function return the exons in
    reversed order. In this case, the first exon relative the the direction of
    transcription (which is probably what the user means, is the last exon
    along the chromosome reading from left to right along the '+' strand (which
    is how the data are stored in UCSC tables).

    E.g.:

        >>> exons({'exonStarts': '10,15', 'exonEnds': '20,25', 'strand': '+'})
        [(10, 20), (15, 25)]
        >>> exons({'exonStarts': '1,3', 'exonEnds': '2,4', 'strand': '-'})
        [(3, 4), (1, 2)]

    Raises a FormattingError when the `row` does not appear to come from a
    valid UCSC gene table.

    """
    try:
        exon_starts = row['exonStarts'].split(',')
        exon_ends = row['exonEnds'].split(',')
    except KeyError as error:
        raise FormattingError(
                "key {!s} not in fields: {!r}".format(
                    error, list(row.keys())))
    positions = []
    for start, end in zip(exon_starts, exon_ends):
        if start != '' and end != '':
            try:
                start, end = int(start), int(end)
            except ValueError:
                raise FormattingError(
                        "unexpected values for 'exonStarts' and 'exonEnds' "
                        "fields: {!r} and {!r}".format(
                            row['exonStarts'], row['exonEnds']))
            if start >= end:
                raise FormattingError(
                        "unexpected values for 'exonStarts' and 'exonEnds' "
                        "fields: {!r} and {!r}. exonEnds values must be "
                        "strictly greater than exonStarts".format(
                            row['exonStarts'], row['exonEnds']))
            positions.append((start, end))
    if row['strand'] == '-':
        positions.reverse()
    return positions


class FormattingError(Exception):
    """Raised when a UCSC file is improperly formatted.

    """
