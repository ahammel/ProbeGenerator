import unittest
import os

from probe_generator import annotation
from probe_generator.sequence import SequenceRange
from probe_generator.test.test_constants import VALIDATION_DATA_DIR

MOCK_ANNOTATION_FILE = [ # input is any iterable of strings
        # UCSC annotation files have a header in this format:
        "#foo\tbar\tbaz\n",
        "ding\tdong\tdong\n",
        ]

MOCK_ROW = {'foo': 'ding', 'bar': 'dong', 'baz': 'dong'}

# Validation test data consists of tables downloaded from UCSC with a mock gene
# 'spiked in'
MOCK_REFSEQ_GENES_FILE = os.path.join(
        VALIDATION_DATA_DIR, "test_refseq_genes.txt")
MOCK_UCSC_GENES_FILE = os.path.join(
        VALIDATION_DATA_DIR, "test_ucsc_genes.txt")


class TestAnnotationValidation(unittest.TestCase):
    """Validation tests for the `annotation` module.

    Validation consists of attempting to retrieve the exons of a mock gene that
    has been manually 'spiked in' to production data. The gene name will always
    be 'MOCK_GENE', and there will always be exactly one exon from bases 2 to 3.

    """
    def assert_mock_gene_in_file(self, annotation_file):
        """Assert that the mock gene is found in the `annotation_file`.

        `annotation_file` is a handle to one of the UCSC annotation files used
        for testing.

        """
        annotations = annotation.parse_ucsc_file(annotation_file)
        matching_features = annotation.lookup_gene("MOCK_GENE", annotations)
        try:
            mock_row, = tuple(matching_features)
            exons = mock_row.exons()
            self.assertEqual(exons, [SequenceRange('0', 2, 3)])
        except ValueError as error:
            self.fail("Unexpected number of mock genes: {}".format(error))

    @unittest.skipIf(not os.path.exists(MOCK_REFSEQ_GENES_FILE),
                     "{} not reachable".format(MOCK_REFSEQ_GENES_FILE))
    def test_refseq_gene_table(self):
        with open(MOCK_REFSEQ_GENES_FILE) as handle:
            self.assert_mock_gene_in_file(handle)

    @unittest.skipIf(not os.path.exists(MOCK_UCSC_GENES_FILE),
                     "{} not reachable".format(MOCK_UCSC_GENES_FILE))
    def test_ucsc_gene_table(self):
        with open(MOCK_UCSC_GENES_FILE) as handle:
            self.assert_mock_gene_in_file(handle)
