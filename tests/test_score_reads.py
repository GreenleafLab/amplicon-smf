import unittest
import numpy as np
from dSMF_footprints_clustering_py3 import score_read_against_region


class TestScoreReadAgainstRegion(unittest.TestCase):
    def test_no_overlap_returns_none(self):
        """If the read does not overlap the region, return (None, None)."""
        L = 10
        region_start = 100
        region_len = L

        rstart = region_start + region_len + 5  # clearly after the region
        rseq = "ACGTACGT"

        mask_CG = np.zeros(L - 1, dtype=bool)
        mask_GC = np.zeros(L - 1, dtype=bool)
        mask_GCG = np.zeros(max(L - 2, 0), dtype=bool)
        mask_C = np.zeros(L, dtype=bool)

        scores, covered = score_read_against_region(
            rstart,
            rseq,
            region_start,
            region_len,
            c_type="CG",
            no_endog=True,
            mask_CG=mask_CG,
            mask_GC=mask_GC,
            mask_GCG=mask_GCG,
            mask_C=mask_C,
        )

        self.assertIsNone(scores)
        self.assertIsNone(covered)

    def test_simple_CG_match(self):
        """CG context: matching CG calls 0 at both positions."""
        ref = "ACGTA"              # CG at positions 1–2
        L = len(ref)
        region_start = 100
        region_len = L
        rstart = region_start
        rseq = ref  # perfect match

        # CG at index 1
        mask_CG = np.array([False, True, False, False], dtype=bool)
        mask_GC = np.zeros(L - 1, dtype=bool)
        mask_GCG = np.zeros(max(L - 2, 0), dtype=bool)
        mask_C = np.array([b == "C" for b in ref], dtype=bool)

        scores, covered = score_read_against_region(
            rstart,
            rseq,
            region_start,
            region_len,
            c_type="CG",
            no_endog=True,
            mask_CG=mask_CG,
            mask_GC=mask_GC,
            mask_GCG=mask_GCG,
            mask_C=mask_C,
        )

        self.assertIsNotNone(scores)
        self.assertEqual(scores.shape[0], L)
        # All covered
        self.assertTrue(np.all(covered))

        # Expect 0 at positions 1 and 2, -1 elsewhere
        expected = np.full(L, -1, dtype=np.int8)
        expected[1] = 0
        expected[2] = 0
        np.testing.assert_array_equal(scores, expected)

    def test_simple_CG_mismatch(self):
        """CG context: non-CG in read at CG site calls 1."""
        ref = "ACGTA"              # CG at positions 1–2
        L = len(ref)
        region_start = 100
        region_len = L
        rstart = region_start

        # mutate C→T so dinucleotide is TG instead of CG
        rseq = "ATGTA"

        mask_CG = np.array([False, True, False, False], dtype=bool)
        mask_GC = np.zeros(L - 1, dtype=bool)
        mask_GCG = np.zeros(max(L - 2, 0), dtype=bool)
        mask_C = np.array([b == "C" for b in ref], dtype=bool)

        scores, covered = score_read_against_region(
            rstart,
            rseq,
            region_start,
            region_len,
            c_type="CG",
            no_endog=True,
            mask_CG=mask_CG,
            mask_GC=mask_GC,
            mask_GCG=mask_GCG,
            mask_C=mask_C,
        )

        expected = np.full(L, -1, dtype=np.int8)
        expected[1] = 1
        expected[2] = 1
        np.testing.assert_array_equal(scores, expected)

    def test_GC_skips_GCG_when_no_endog_false(self):
        """
        GC context: when no_endog is False, positions that are part of GCG
        (marked in mask_GCG) are skipped.
        """
        # ref indices: 0 1 2
        #             G C G
        # GC at (0,1) and (1,2), but both are part of GCG starting at 0
        ref = "GCG"
        L = len(ref)
        region_start = 100
        region_len = L
        rstart = region_start
        rseq = ref  # matching

        # GC starts at 0 and 1 for this short ref
        mask_GC = np.array([True, True], dtype=bool)
        # GCG starts at index 0
        mask_GCG = np.array([True], dtype=bool)   # same indexing convention for start
        mask_CG = np.zeros(L - 1, dtype=bool)
        mask_C = np.array([b == "C" for b in ref], dtype=bool)

        scores, covered = score_read_against_region(
            rstart,
            rseq,
            region_start,
            region_len,
            c_type="GC",
            no_endog=False,        # => skip_gcg=True
            mask_CG=mask_CG,
            mask_GC=mask_GC,
            mask_GCG=mask_GCG,
            mask_C=mask_C,
        )

        # Because both GC positions are part of GCG, they should be skipped,
        # leaving all scores = -1.
        expected = np.full(L, -1, dtype=np.int8)
        np.testing.assert_array_equal(scores, expected)

    def test_GC_no_skip_when_no_endog_true(self):
        """
        GC context: when no_endog is True, GCG positions are NOT skipped.
        """
        ref = "GCG"
        L = len(ref)
        region_start = 100
        region_len = L
        rstart = region_start
        rseq = ref

        mask_GC = np.array([True, True], dtype=bool)
        mask_GCG = np.array([True], dtype=bool)
        mask_CG = np.zeros(L - 1, dtype=bool)
        mask_C = np.array([b == "C" for b in ref], dtype=bool)

        scores, covered = score_read_against_region(
            rstart,
            rseq,
            region_start,
            region_len,
            c_type="GC",
            no_endog=True,         # => skip_gcg=False
            mask_CG=mask_CG,
            mask_GC=mask_GC,
            mask_GCG=mask_GCG,
            mask_C=mask_C,
        )

        # Two GC start positions; each GC scores 0 at i and i+1
        expected = np.full(L, -1, dtype=np.int8)
        # GC at 0–1 -> indices 0,1
        expected[0] = 0
        expected[1] = 0
        # GC at 1–2 -> indices 1,2
        expected[1] = 0
        expected[2] = 0
        np.testing.assert_array_equal(scores, expected)

    def test_allC_scoring(self):
        """allC: every C gets scored 0 if read has C, 1 otherwise."""
        ref = "ACCA"               # Cs at positions 1 and 2
        L = len(ref)
        region_start = 100
        region_len = L
        rstart = region_start

        # Case 1: read matches ref => all Cs are unconverted (0)
        rseq1 = ref
        mask_C = np.array([b == "C" for b in ref], dtype=bool)
        mask_CG = np.zeros(L - 1, dtype=bool)
        mask_GC = np.zeros(L - 1, dtype=bool)
        mask_GCG = np.zeros(max(L - 2, 0), dtype=bool)

        scores1, covered1 = score_read_against_region(
            rstart,
            rseq1,
            region_start,
            region_len,
            c_type="allC",
            no_endog=True,
            mask_CG=mask_CG,
            mask_GC=mask_GC,
            mask_GCG=mask_GCG,
            mask_C=mask_C,
        )

        expected1 = np.full(L, -1, dtype=np.int8)
        expected1[1] = 0
        expected1[2] = 0
        np.testing.assert_array_equal(scores1, expected1)

        # Case 2: one C converted to T => that position scores 1
        rseq2 = "ATCA"  # positions: A T C A; C only at index 2
        scores2, covered2 = score_read_against_region(
            rstart,
            rseq2,
            region_start,
            region_len,
            c_type="allC",
            no_endog=True,
            mask_CG=mask_CG,
            mask_GC=mask_GC,
            mask_GCG=mask_GCG,
            mask_C=mask_C,
        )

        expected2 = np.full(L, -1, dtype=np.int8)
        # index 1: ref C, read T -> 1
        expected2[1] = 1
        # index 2: ref C, read C -> 0
        expected2[2] = 0
        np.testing.assert_array_equal(scores2, expected2)


if __name__ == "__main__":
    unittest.main()