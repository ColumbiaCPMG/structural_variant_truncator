import unittest
import sv_truncator as T
import portion as P


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

        Result (1 SV spanning first and last exons):
                =======================================
                ^100                                  ^300
        """
        regions = {}
        regions["chr1"] = [
            P.closed(100, 200),
            P.closed(150, 250),
            P.closed(275, 400),
        ]
        sv = P.closed(50, 300)

        result = T.getMatchingIntervals(
            "chr1", sv, regions
        )

        intervals = P.closed(100, 250) | P.closed(275, 300)
        self.assertEqual(
            result,
            intervals
        )

    def test_chromosome_does_not_exist_in_regions_file(self):
        regions = {}
        regions["chr1"] = ''
        sv = P.closed(0, 1)

        result = T.getMatchingIntervals(
            "chr10", sv, regions
        )
        self.assertEqual(result, P.empty())

    def test_chromosome_already_prefixed_with_chr_in_regions_file(self):
        regions, trees = T.loadGenomicCoordinatesFile(
                'test_files/Homo_sapiens.GRCh37.75_10_chr.txt')

        sv = P.closed(10, 76)
        result = T.getMatchingIntervals(
            "chr1", sv, regions
        )

        expected_interval = P.closed(10, 45) | P.closed(47, 75)
        self.assertEqual(result, expected_interval)

    def test_positions_are_not_numeric_in_regions_file(self):
        pass

    def test_positions_are_not_numeric_in_PennCNV_file(self):
        pass

    def test_input_file_has_blank_lines(self):
        pass






if __name__ == '__main__':
    unittest.main()
