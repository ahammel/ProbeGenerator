import unittest
import sys
import io

import mock

from probe_generator import sequence, probe_statement


def mock_exons(feature):
    """Mock exon function.

    """
    if feature['gene'] == 'FOO':
        return [(50, 150)]
    elif feature['gene'] == 'BAR':
        return [(200, 300), (20, 62)]
    elif feature['gene'] == 'LEFT':
        return [(100, 200)]
    elif feature['gene'] == 'RIGHT':
        return[(300, 400)]


class AbstractSequenceRangeTestCase(unittest.TestCase):
    """Provide setUp function for `sequence.sequence_range` tests.

    """
    def setUp(self):
        self.patcher = mock.patch(
                'probe_generator.annotation.exons',
                mock_exons)
        self.patcher.start()

    def tearDown(self):
        self.patcher.stop()


class TestSequenceRangePositional(AbstractSequenceRangeTestCase):
    """Test cases for the sequence.sequence_range method using the
    `positional` strategy.


    """
    def setUp(self):
        super().setUp()
        self.probe_specification = {
                'gene1': 'FOO',
                'feature1': ('exon', 1),
                'side1': 'start',
                'bases1': 50,
                'gene2': 'BAR',
                'feature2': ('exon', 2),
                'side2': 'end',
                'bases2': 20,
                'separator': '/',
                }
        self.feature_1 = {
                'gene': 'FOO', # Fake field for mocking purposes
                'strand': '+',
                'chrom': 'chr1',
                'exonStarts': '50,',
                'exonEnds':  '150,'
                }
        self.feature_2 = {
                'gene': 'BAR',
                'strand': '-',
                'chrom': 'chr2',
                'exonStarts': '20,200,', # Feature is on '-' strand, so exon
                'exonEnds':  '62,300'    # numbers are reversed relative to '+'
                }

    def test_sequence_range_returns_incusive_ranges_and_rev_comp_flag(self):
        """
        The sequence_range function returns a list of pairs of tuples
        specifying the range of base pairs to be extracted from the genome
        with a flag specifying whether the second sequence should be
        reverse-complemented.

        """
        range_spec = sequence.sequence_range(
                self.probe_specification,
                self.feature_1,
                self.feature_2)
        self.assertEqual(
                range_spec,
                {'chromosome1': '1',
                 'start1':      50,
                 'end1':        99,
                 'chromosome2': '2',
                 'start2':      20,  # Feature is on '-' strand, so
                 'end2':        39,  # start and end are reversed
                 'inversion':   True})

    def test_sequence_range_returns_one_base_range_when_bases_is_1(self):
        self.probe_specification['bases1'] = 1
        self.probe_specification['bases2'] = 1
        # Both modified for inherited tests
        range_spec = sequence.sequence_range(
                self.probe_specification,
                self.feature_1,
                self.feature_2)
        self.assertEqual(
                (range_spec['start1'], range_spec['end1']),
                (50, 50))

    def test_sequence_range_returns_entire_feature_when_bases_globbed(self):
        self.probe_specification['bases1'] = '*'
        self.probe_specification['bases2'] = '*'
        # Both globbed for inherited tests
        range_spec = sequence.sequence_range(
                self.probe_specification,
                self.feature_1,
                self.feature_2)
        self.assertEqual(
                (range_spec['chromosome1'],
                 range_spec['start1'],
                 range_spec['end1']),
                ('1', 50, 150))

    def test_raises_InterfaceError_when_specification_missing_feature(self):
        with self.assertRaises(sequence.InterfaceError):
            del self.probe_specification['feature1']
            sequence.sequence_range(
                    self.probe_specification,
                    self.feature_1,
                    self.feature_2)

    def test_raises_InterfaceError_when_specification_missing_bases(self):
        with self.assertRaises(sequence.InterfaceError):
            del self.probe_specification['bases1']
            sequence.sequence_range(
                    self.probe_specification,
                    self.feature_1,
                    self.feature_2)

    def test_raises_InterfaceError_when_specification_missing_side(self):
        with self.assertRaises(sequence.InterfaceError):
            del self.probe_specification['side1']
            sequence.sequence_range(
                    self.probe_specification,
                    self.feature_1,
                    self.feature_2)

    def test_raises_InterfaceError_when_specification_missing_separator(self):
        with self.assertRaises(sequence.InterfaceError):
            del self.probe_specification['separator']
            sequence.sequence_range(
                    self.probe_specification,
                    self.feature_1,
                    self.feature_2)

    def test_raises_InterfaceError_when_row_missing_strand(self):
        with self.assertRaises(sequence.InterfaceError):
            del self.feature_1['strand']
            print(sequence.sequence_range(
                    self.probe_specification,
                    self.feature_1,
                    self.feature_2))

    def test_raises_NoFeatureError_when_exon_out_of_range(self):
        message = ("specification requires feature 'exon'\[20\], but "
                   "row specifies only [12] 'exon'\(s\)")
        # Message matches for the inherited tests in
        with self.assertRaisesRegex(sequence.NoFeatureError, message):
            self.probe_specification['feature1'] = ('exon', 20)
            sequence.sequence_range(
                    self.probe_specification,
                    self.feature_1,
                    self.feature_2)


class TestSequenceRangesEndwise(AbstractSequenceRangeTestCase):
    """Exhaustive tests for all combinations of sides and strands for
    `sequence.sequence_range` using the positional strategy.

    """
    def setUp(self):
        super().setUp()
        self.endwise_specification = {
                'gene1': 'LEFT',
                'feature1': ('exon', 1),
                'side1': None, # must be set in test
                'bases1': 50,
                'gene2': 'RIGHT',
                'feature2': ('exon', 1),
                'side2': None,
                'bases2': 50,
                'separator': '/',
                }
        self.endwise_left = {
                'gene': 'LEFT',
                'strand': None, # set in test
                'chrom': 'chr1',
                'exonStarts': '100,',
                'exonEnds':  '200,'
                }
        self.endwise_right = {
                'gene': 'RIGHT',
                'strand': None, # set in test
                'chrom': 'chr2',
                'exonStarts': '300,',
                'exonEnds':  '400,'
                }

    def endwise_test(self, sides_and_strands, result):
        """Assert that the results are correct, given the sides and strands of
        the features.

        `sides_and_strands` is the value of the side of the left feature, the
        side of the right feature, the strand of the left feature and the
        strand of the right feature in that order.

        These values are assigned to the `endwise_specification`,
        `endwise_left` and `endwise_right`, which are then passed to
        `sequence_ranges`.

        The return values is compared to `result`, which has the chromosomes
        specified in the test for ease of typing.

        """
        left_side, right_side, left_strand, right_strand = sides_and_strands
        self.endwise_specification['side1'] = left_side
        self.endwise_specification['side2'] = right_side
        self.endwise_left['strand'] = left_strand
        self.endwise_right['strand'] = right_strand
        the_result = dict(result, chromosome1='1', chromosome2='2')
        self.assertEqual(
                sequence.sequence_range(
                    self.endwise_specification,
                    self.endwise_left,
                    self.endwise_right),
                the_result)

    def test_endwise_start_start_plus_plus(self):
        self.endwise_test(
                ('start', 'start', '+', '+'),
                {'start1': 100, 'end1': 149,
                 'start2': 300, 'end2': 349,
                 'inversion': True})

    def test_endwise_start_start_plus_minus(self):
        self.endwise_test(
                ('start', 'start', '+', '-'),
                {'start1': 100, 'end1': 149,
                 'start2': 351, 'end2': 400,
                 'inversion': False})

    def test_endwise_start_start_minus_plus(self):
        self.endwise_test(
                ('start', 'start', '-', '+'),
                {'start1': 151, 'end1': 200,
                 'start2': 300, 'end2': 349,
                 'inversion': False})

    def test_endwise_start_start_minus_minus(self):
        self.endwise_test(
                ('start', 'start', '-', '-'),
                {'start1': 151, 'end1': 200,
                 'start2': 351, 'end2': 400,
                 'inversion': True})

    def test_endwise_start_end_plus_plus(self):
        self.endwise_test(
                ('start', 'end', '+', '+'),
                {'start1': 100, 'end1': 149,
                 'start2': 351, 'end2': 400,
                 'inversion': False})

    def test_endwise_start_end_plus_minus(self):
        self.endwise_test(
                ('start', 'end', '+', '-'),
                {'start1': 100, 'end1': 149,
                 'start2': 300, 'end2': 349,
                 'inversion': True})

    def test_endwise_start_end_minus_plus(self):
        self.endwise_test(
                ('start', 'end', '-', '+'),
                {'start1': 151, 'end1': 200,
                 'start2': 351, 'end2': 400,
                 'inversion': True})

    def test_endwise_start_end_minus_minus(self):
        self.endwise_test(
                ('start', 'end', '-', '-'),
                {'start1': 151, 'end1': 200,
                 'start2': 300, 'end2': 349,
                 'inversion': False})

    def test_endwise_end_start_plus_plus(self):
        self.endwise_test(
                ('end', 'start', '+', '+'),
                {'start1': 151, 'end1': 200,
                 'start2': 300, 'end2': 349,
                 'inversion': False})

    def test_endwise_end_start_plus_minus(self):
        self.endwise_test(
                ('end', 'start', '+', '-'),
                {'start1': 151, 'end1': 200,
                 'start2': 351, 'end2': 400,
                 'inversion': True})

    def test_endwise_end_start_minus_plus(self):
        self.endwise_test(
                ('end', 'start', '-', '+'),
                {'start1': 100, 'end1': 149,
                 'start2': 300, 'end2': 349,
                 'inversion': True})

    def test_endwise_end_start_minus_minus(self):
        self.endwise_test(
                ('end', 'start', '-', '-'),
                {'start1': 100, 'end1': 149,
                 'start2': 351, 'end2': 400,
                 'inversion': False})

    def test_endwise_end_end_plus_plus(self):
        self.endwise_test(
                ('end', 'end', '+', '+'),
                {'start1': 151, 'end1': 200,
                 'start2': 351, 'end2': 400,
                 'inversion': True})

    def test_endwise_end_end_plus_minus(self):
        self.endwise_test(
                ('end', 'end', '+', '-'),
                {'start1': 151, 'end1': 200,
                 'start2': 300, 'end2': 349,
                 'inversion': False})

    def test_endwise_end_end_minus_plus(self):
        self.endwise_test(
                ('end', 'end', '-', '+'),
                {'start1': 100, 'end1': 149,
                 'start2': 351, 'end2': 400,
                 'inversion': False})

    def test_endwise_end_end_minus_minus(self):
        self.endwise_test(
                ('end', 'end', '-', '-'),
                {'start1': 100, 'end1': 149,
                 'start2': 300, 'end2': 349,
                 'inversion': True})


class TestSequenceRangeReadThrough(TestSequenceRangePositional):
    """Test cases for the read-through strategy for sequence.sequence_range

    Runs the tests from the TestSequenceRangePositional suite, but with a
    suitable read-through specification substituted.

    There are additional tests for the various combinations of strands, and
    for the warnings that should be raised when the sides of the specification
    don't make sense.


    """
    def setUp(self):
        super().setUp()
        self.probe_specification = {
                'gene1':     'BAR',
                'feature1':  ('exon', 2),
                'side1':     'end',
                'bases1':    20,
                'gene2':     'FOO',
                'feature2':  ('exon', 1),
                'side2':     'start',
                'bases2':    50,
                'separator': '->',
                }
        self.feature_1, self.feature_2 = self.feature_2, self.feature_1

        self.stderr_backup = sys.stderr
        sys.stderr = io.StringIO()

    def tearDown(self):
        super().tearDown()
        sys.stderr.close()
        sys.stderr = self.stderr_backup

    def chromosomes_swapped(self):
        """Return True if the chromosome1 of the coordinate specification
        returned by sequence.sequence_range corresponds to
        self.feature2['chrom'] and vice versa.

        """
        coords = sequence.sequence_range(
                    self.probe_specification,
                    self.feature_1,
                    self.feature_2)
        chromosome1, chromosome2 = (self.feature_1['chrom'],
                                    self.feature_2['chrom'])
        return ((coords['chromosome1'], coords['chromosome2']) ==
                (chromosome2.lstrip('chr'), chromosome1.lstrip('chr')))

    def assert_warning_raised(self, warning):
        """Assert that the `warning` was printed to stderr when
        sequence.sequence_range was called with the standard parameters

        """
        sequence.sequence_range(
                self.probe_specification,
                self.feature_1,
                self.feature_2)
        message = sys.stderr.getvalue()
        self.assertEqual(message, warning)

    def test_genes_not_swapped_when_both_features_on_plus_strand(self):
        self.feature_1['strand'] = '+'
        self.feature_2['strand'] = '+'
        self.assertFalse(self.chromosomes_swapped())

    def test_genes_not_swapped_when_feature_1_plus_feature_2_minus(self):
        self.feature_1['strand'] = '+'
        self.feature_2['strand'] = '-'
        self.assertFalse(self.chromosomes_swapped())

    def test_genes_swapped_when_feature_1_minus_feature_2_plus(self):
        self.feature_1['strand'] = '-'
        self.feature_2['strand'] = '+'
        self.assertTrue(self.chromosomes_swapped())

    def test_genes_swapped_when_feature_1_minus_feature_2_minus(self):
        self.feature_1['strand'] = '-'
        self.feature_2['strand'] = '-'
        self.assertTrue(self.chromosomes_swapped())

    def test_no_warning_raised_when_side_1_is_end_and_side_2_is_start(self):
        # this is the case by default
        self.assert_warning_raised('')

    def test_warning_raised_when_side_1_is_start(self):
        self.probe_specification['side1'] = 'start'
        self.assert_warning_raised(sequence.WARNING_MESSAGE)

    def test_warning_raised_when_side_2_is_end(self):
        self.probe_specification['side2'] = 'end'
        self.assert_warning_raised(sequence.WARNING_MESSAGE)

    def test_warning_raised_when_side_1_is_staat_and_side_2_is_end(self):
        self.probe_specification['side1'] = 'start'
        self.probe_specification['side2'] = 'end'
        self.assert_warning_raised(sequence.WARNING_MESSAGE)


class TestExonsIntegrationPositional(TestSequenceRangePositional):
    """Integration tests for `sequence.sequence_range` and `annotation.exons`
    with a positional specification.

    Same tests as TestSequenceRanges with the production `annotation.exons`
    function instead of the mock.

    """
    def setUp(self):
        super().setUp()
        mock.patch.stopall()

    def tearDown(self):
        pass


class TestExonsIntegrationReadThrough(TestSequenceRangeReadThrough):
    """Integration tests for `sequence.sequence_range` and `annotation.exons`
    with a read-through specification.

    Same tests as TestSequenceRanges with the production `annotation.exons`
    function instead of the mock.

    """
    def setUp(self):
        super().setUp()
        mock.patch.stopall()

    def tearDown(self):
        sys.stderr.close()
        sys.stderr = self.stderr_backup


class TestProbeStatementParseIntegrationPositional(TestSequenceRangePositional):
    """Integration tests for `sequence.sequence_range` and `probe_statement.parse`
    with a positional specification.

    Same tests as `TestSequenceRanges`, but with the `probe_specification`
    parsed from a probe statement instead of mocked in.

    """
    def setUp(self):
        super().setUp()
        self.probe_specification = probe_statement.parse(
                "FOO#exon[1] +50 / BAR#exon[2] -20")


class TestProbeStatementParseIntegrationReadThrough(TestSequenceRangeReadThrough):
    """Integration tests for `sequence.sequence_range` and `probe_statement.parse`
    with a read-through specification.

    Same tests as `TestSequenceRanges`, but with the `probe_specification`
    parsed from a probe statement instead of mocked in.

    """
    def setUp(self):
        super().setUp()
        self.probe_specification = probe_statement.parse(
                "BAR#exon[2] -20 -> FOO#exon[1] + 50")
