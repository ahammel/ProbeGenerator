import unittest

import mock

from probe_generator import sequence, probe_statement


def mock_exons(feature):
    """Mock exon function.

    """
    if feature['feature_number_sentinel'] == 1:
        return [(50, 150)]
    elif feature['feature_number_sentinel'] == 2:
        return [(1, 10), (20, 62)]


class TestSequenceRanges(unittest.TestCase):
    """Test cases for the sequence.sequence_range method.

    """
    def setUp(self):
        self.probe_specification = {
                'gene1': 'FOO',
                'gene2': 'BAR',
                'feature1': ('exon', 1),
                'feature2': ('exon', 2),
                'side1': 'start',
                'side2': 'end',
                'bases1': 50,
                'bases2': 20,
                }
        self.feature_1 = {
                'feature_number_sentinel': 1, # Fake field for mocking purposes
                'strand': '+',
                'chromosome': 'chr1',
                'exonStarts': '50,',
                'exonEnds':  '150,'
                }
        self.feature_2 = {
                'feature_number_sentinel': 2,
                'strand': '-',
                'chromosome': 'chr2',
                'exonStarts': '1,20,',
                'exonEnds':  '10,62,'
                }

        self.patcher = mock.patch(
                'probe_generator.annotation.exons',
                mock_exons)
        self.patcher.start()

    def tearDown(self):
        self.patcher.stop()

    def assert_reverse_complement_flag_equals(self, flag):
        """Call `sequence.sequence_range` on the internal feature and
        probe_specification objects and assert that the 'reverse' value of the
        return value is equal to 'flag'.

        """
        range_spec = sequence.sequence_range(
                self.probe_specification,
                self.feature_1,
                self.feature_2)
        self.assertEqual(range_spec['reverse'], flag)


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
                range_spec['range_1'],
                ('1', 50, 100))
        self.assertEqual(
                range_spec['range_2'],
                ('2', 42, 62))
        self.assertTrue('reverse' in range_spec)

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
