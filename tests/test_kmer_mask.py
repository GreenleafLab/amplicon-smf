import unittest
import numpy as np
from numpy.testing import assert_array_equal
from dSMF_footprints_clustering_py3 import kmer_mask  # adjust if needed


class TestKmerMask(unittest.TestCase):
    def test_single_hit_monomer(self):
        seq = "TGCGAC"
        kmer = "A"
        out = kmer_mask(seq, kmer)
        assert_array_equal(out, np.array([False, False, False, False, True, False], dtype=bool))
        self.assertEqual(len(out), len(seq) - len(kmer) + 1)
        self.assertEqual(out.dtype, bool)

    def test_single_hit_dimer(self):
        out = kmer_mask("ATGC", "TG")
        assert_array_equal(out, np.array([False, True, False], dtype=bool))
        self.assertEqual(len(out), len("ATGC") - len("TG") + 1)
        self.assertEqual(out.dtype, bool)

    def test_overlapping_hits(self):
        # "AAAAA" with "AAA" hits at positions 0,1,2
        out = kmer_mask("AAAAA", "AAA")
        assert_array_equal(out, np.array([True, True, True], dtype=bool))

    def test_no_hits(self):
        out = kmer_mask("ACGTACGT", "TTT")
        assert_array_equal(out, np.zeros(len("ACGTACGT") - 3 + 1, dtype=bool))

    def test_seq_shorter_than_kmer(self):
        out = kmer_mask("AC", "ACG")
        self.assertEqual(out.size, 0)
        self.assertEqual(out.dtype, bool)

    def test_length_formula_parametric(self):
        cases = [
            ("A", "A", 1),
            ("AC", "C", 2),
            ("ACG", "ACG", 1),
            ("ACGT", "CG", 3),
            ("TTTT", "TT", 3),
        ]
        for seq, kmer, _ in cases:
            with self.subTest(seq=seq, kmer=kmer):
                out = kmer_mask(seq, kmer)
                expected_len = max(0, len(seq) - len(kmer) + 1)
                self.assertEqual(len(out), expected_len)

    def test_empty_zeromer(self):
        # Current implementation is all True.
        out = kmer_mask("ACG", "")
        assert_array_equal(out, np.ones(len("ACG") + 1, dtype=bool))

    def test_non_ascii_raises(self):
        with self.assertRaises(UnicodeEncodeError):
            kmer_mask("AÄC", "A")  # function encodes to ASCII


if __name__ == "__main__":
    unittest.main()
