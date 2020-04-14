import unittest
import sv_truncator as T
from intervaltree import Interval, IntervalTree


class Tests(unittest.TestCase):
    def test_sv_truncated_with_three_exons(self):
        """
        Test one large SV overlapping 3 exons truncated to first and last exon

        Large SV:
        ===============================================
        ^50                                           ^300

        Smaller exonic regions (100->200, 150->250, 275->400):
                ====================
                      ======================
                                                  =====================
                ^100                                                  ^400

        Result (1 SV spanning first and last exon):
                =======================================
                ^100                                  ^300
        """
        regions = {}
        regions["chr1"] = IntervalTree(
            [
                Interval(100, 200),
                Interval(150, 250),
                Interval(275, 400)
            ]
        )

        sv = Interval(50, 300)

        result = T.getMatchingIntervalsFromTree(
            "chr1", sv, regions
        )
        expected = IntervalTree([Interval(100, 300)])

        self.assertEqual(
            result,
            expected
        )

    def test_chromosome_does_not_exist_in_regions_file(self):
        trees = {}
        trees['chr1'] = IntervalTree([Interval(-99, 99)])
        sv = Interval(10, 45)

        result = T.getMatchingIntervalsFromTree(
            "chr10", sv, trees
        )
        self.assertEqual(result.is_empty(), IntervalTree().is_empty())

    def test_chromosome_already_prefixed_with_chr_in_regions_file(self):
        trees = T.loadGenomicCoordinatesFile(
            'test_files/Homo_sapiens.GRCh37.75_10_chr.txt')

        sv = Interval(12, 76)
        result = T.getMatchingIntervalsFromTree(
            "chr1", sv, trees
        )

        expected_interval = IntervalTree([Interval(12, 76)])
        self.assertEqual(result, expected_interval)

    def test_no_matching_region_returns_empty_tree(self):
        trees = {}
        trees['chr1'] = IntervalTree([Interval(0, 10)])
        sv = Interval(12, 76)
        result = T.getMatchingIntervalsFromTree(
            "chr1", sv, trees
        )

        expected = IntervalTree()
        self.assertEqual(result, expected)


    # def test_input_file_has_blank_lines(self):
    #     pass

    # def test_positions_are_not_numeric_in_regions_file(self):
    #     pass

    # def test_positions_are_not_numeric_in_PennCNV_file(self):
    #     pass


if __name__ == '__main__':
    unittest.main()
