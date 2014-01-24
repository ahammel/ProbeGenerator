import unittest

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


class TestSequenceRanges(unittest.TestCase):
    """Test cases for the sequence.sequence_range method.

    """
    def setUp(self):
        self.probe_specification = {
                'gene1': 'FOO',
                'feature1': ('exon', 1),
                'side1': 'start',
                'bases1': 50,
                'gene2': 'BAR',
                'feature2': ('exon', 2),
                'side2': 'end',
                'bases2': 20,
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

        # The specification and features that follow are for the 'endwise'
        # family of tests, which exhaustively test sequence_ranges under all
        # possible combinations of orientations and feature ends.
        self.endwise_specification = {
                'gene1': 'LEFT',
                'feature1': ('exon', 1),
                'side1': None, # must be set in test
                'bases1': 50,
                'gene2': 'RIGHT',
                'feature2': ('exon', 1),
                'side2': None,
                'bases2': 50,
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

        # `sequence_ranges` calls `exons`
        self.patcher = mock.patch(
                'probe_generator.annotation.exons',
                mock_exons)
        self.patcher.start()

    def tearDown(self):
        self.patcher.stop()

    def assert_reverse_complement_flag_equals(self, flag):
        """Call `sequence.sequence_range` on the internal feature and
        probe_specification objects and assert that the 'inversion' value of the
        return value is equal to 'flag'.

        """
        range_spec = sequence.sequence_range(
                self.probe_specification,
                self.feature_1,
                self.feature_2)
        self.assertEqual(range_spec['inversion'], flag)

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
        range_spec = sequence.sequence_range(
                self.probe_specification,
                self.feature_1,
                self.feature_2)
        self.assertEqual(
                (range_spec['start1'], range_spec['end1']),
                (50, 50))

    def test_sequence_range_returns_entire_feature_when_bases_globbed(self):
        self.probe_specification['bases1'] = '*'
        range_spec = sequence.sequence_range(
                self.probe_specification,
                self.feature_1,
                self.feature_2)
        self.assertEqual(
                (range_spec['chromosome1'],
                 range_spec['start1'],
                 range_spec['end1']),
                ('1', 50, 150))

    def test_rev_comp_True_when_strands_opposite_and_sides_opposite(self):
        # This is the case by default
        self.assert_reverse_complement_flag_equals(True)

    def test_rev_comp_False_when_strands_equal_and_sides_opposite(self):
        self.feature_1['strand'] = '-'
        self.assert_reverse_complement_flag_equals(False)

    def test_rev_comp_False_when_strands_opposite_and_sides_equal(self):
        self.probe_specification['side1'] = 'end'
        self.assert_reverse_complement_flag_equals(False)

    def test_rev_comp_True_when_strands_equal_and_sides_equal(self):
        self.feature_1['strand'] = '-'
        self.probe_specification['side1'] = 'end'
        self.assert_reverse_complement_flag_equals(True)

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

    def test_raises_NoFeatureError_when_exon_out_of_range(self):
        message = ("specification requires feature 'exon'\[20\], but "
                   "row specifies only 1 'exon'\(s\)")
        with self.assertRaisesRegex(sequence.NoFeatureError, message):
            self.probe_specification['feature1'] = ('exon', 20)
            sequence.sequence_range(
                    self.probe_specification,
                    self.feature_1,
                    self.feature_2)

    def endwise_test(self, sides_and_strands, result):
        """Assert that the results is correct, given the sides and strands of
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


class TestSequenceRangesExonsIntegration(TestSequenceRanges):
    """Integration tests for `sequence.sequence_range` and `annotation.exons`.

    Same tests as TestSequenceRanges with the production `annotation.exons`
    function instead of the mock.

    """
    def setUp(self):
        super(TestSequenceRangesExonsIntegration, self).setUp()
        mock.patch.stopall()

    def tearDown(self):
        pass


class TestSequenceRangesProbeStatementParseIntegration(TestSequenceRanges):
    """Integration tests for `sequence.sequence_range` and `probe_statement.parse`.

    Same tests as `TestSequenceRanges`, but with the `probe_specification`
    parsed from a probe statement instead of mocked in.

    """
    def setUp(self):
        super(TestSequenceRangesProbeStatementParseIntegration, self).setUp()
        self.probe_specification = probe_statement.parse(
                "FOO#exon[1] +50 / BAR#exon[2] -20")
        self.endwise_specification = probe_statement.parse(
                "LEFT#exon[1]*50/RIGHT#exon[1]*50")
