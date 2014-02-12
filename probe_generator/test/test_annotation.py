import unittest
import os
import re

import mock

from probe_generator import annotation
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


class TestAnnotaion(unittest.TestCase):
    """Test cases for annotation parsing.

    """
    def setUp(self):
        self.annotation_file = iter(MOCK_ANNOTATION_FILE)

    def test_mock_annotation_parsed_correctly(self):
        self.assertEqual(
                list(annotation.parse_ucsc_file(MOCK_ANNOTATION_FILE)),
                [MOCK_ROW])


class TestLookupGene(unittest.TestCase):
    """Test cases for getting the data of genes by name.

    """
    def setUp(self):
        self.annotation = iter([MOCK_ROW])

    @mock.patch.object(annotation, "_GENE_NAME_FIELDS", ('foo',))
    def test_lookup_gene_yields_rows_matching_gene_name_for_refseq_file(self):
        """
        By patching the _GENE_NAME_FIELDS module constant, we specify that gene
        names should be looked up under the field named 'foo'.

        """
        self.assertEqual(
                list(annotation.lookup_gene('ding', self.annotation)),
                [MOCK_ROW])

    def test_lookup_gene_returns_empty_generator_for_empty_annotation(self):
        self.assertEqual(
                list(annotation.lookup_gene('wilder', self.annotation)),
                [])


class TestExons(unittest.TestCase):
    """Test cases for the annotation.exons function

    """
    def setUp(self):
        self.row = {
                'strand': '+',
                'exonStarts': '10,30,50,',
                'exonEnds': '20,40,60,'}

    def test_exons_returns_exon_position_tuples(self):
        """
        Given a properly-formatted UCSC table row, annotations.exons returns a
        list of tuples the exon's start position and the exon's end position.

        """
        self.assertEqual(
                annotation.exons(self.row),
                [(10, 20),
                 (30, 40),
                 (50, 60)])

    def test_exons_returns_reversed_positions_when_strand_minus(self):
        self.row['strand'] = '-'
        self.assertEqual(
                annotation.exons(self.row),
                [(50, 60),
                 (30, 40),
                 (10, 20)])

    def test_exons_raises_FormattingError_on_empty_dictionary(self):
        message = "key 'exonStarts' not in fields: \[\]"
        with self.assertRaisesRegex(annotation.FormattingError, message):
            annotation.exons({})

    def test_exons_raises_FormattingError_on_non_int_values_in_string(self):
        message = re.escape(
                   "unexpected values for 'exonStarts' and 'exonEnds' fields: "
                   "'1,banana,' and '2,surprise,'")
        with self.assertRaisesRegex(annotation.FormattingError, message):
            annotation.exons({
                'exonStarts': '1,banana,',
                'exonEnds': '2,surprise,',
                'strand': '-'
                })

    def test_exons_raises_FormattingError_when_start_greater_than_end(self):
        message = re.escape(
                "unexpected values for 'exonStarts' and 'exonEnds' fields: "
                "'10,20,' and '1,2,'. exonEnds values must be strictly "
                "greater than exonStarts")
        with self.assertRaisesRegex(annotation.FormattingError, message):
            annotation.exons({
                'exonStarts': '10,20,',
                'exonEnds': '1,2,',
                'strand': '-'})

    def test_exons_raises_FormattingError_when_row_missing_strand(self):
        message = re.escape(
                "key 'strand' not in fields: ['exonEnds', 'exonStarts']")
        with self.assertRaisesRegex(annotation.FormattingError, message):
            del self.row['strand']
            annotation.exons(self.row)


class TestAnnotationLookupGeneIntegration(TestLookupGene):
    """Integration tests for `annotation.parse_annotation` and
    `annotation.lookup_gene`

    This is done by re-running the TestLookupGene test using the output of
    `parse_annotation` as the fixture instead of a mock.

    """
    def setUp(self):
        self.annotation = annotation.parse_ucsc_file(
                iter(MOCK_ANNOTATION_FILE))


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
            exons = annotation.exons(mock_row)
            self.assertEqual(exons, [(2, 3)])
        except ValueError as error:
            self.fail(
                    "Unexpected number of mock genes: {}".format(error))
        except annotation.FormattingError as error:
            raise annotation.FormattingError(
                    "Mock row formatted improperly: {}".format(error))

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
