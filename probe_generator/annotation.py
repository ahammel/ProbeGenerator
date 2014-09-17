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
    4. Add the gene identifier field to the _GENE_NAME_FIELDS constant in the
       'transcript' module. Make sure the validation test passes.
    5. Add the name of the gene table to list of supported tables above.

"""
import csv

from probe_generator.transcript import Transcript


def parse_ucsc_file(handle):
    """Return a csv.DictReader object relating fields of UCSC annotation files
    to values.

    `handle` is an iterator of the lines of file, which must be in the standard
    UCSC format (i.e., tab-delimited file, first line starts with '#' and
    specifies field names.

    """
    lines = (line.lstrip('#') for line in handle)
    return (Transcript(line) for line in
            csv.DictReader(lines, dialect='excel-tab'))


def lookup_gene(gene_name, ucsc_file):
    """Yield data for features in a `ucsc_file` for a specific gene.

    `ucsc_file` is an iterator of dictionaries giving the data from a UCSC gene
    file, as might be returned by `parse_ucsc_file`. Currently supported
    formats are given in the docstring of the `annotation` module.

    """
    for transcript in ucsc_file:
        if transcript.gene_id == gene_name:
            yield transcript
