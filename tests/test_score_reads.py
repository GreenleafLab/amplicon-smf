import unittest
import numpy as np
from dSMF_footprints_clustering_py3 import score_read_against_region


class TestScoreReadAgainstRegion(unittest.TestCase):
    def test_no_overlap_returns_none(self):
        """Read completely outside region → (None, None)."""
        ref = "ACGTACGT"
        # region is [0, 8); read starts far to the right
        scores, covered = score_read_against_region(
            rstart=100,
            rseq="ACGT",
            region_start=0,
            c_type="CG",
            no_endog=False,
            ref=ref,
        )
        self.assertIsNone(scores)
        self.assertIsNone(covered)

    def test_coverage_mask_partial_overlap(self):
        """Covered positions should match the read overlap with the region."""
        ref = "ACGTACGT"
        # region genomic [100, 108), indices 0..7
        region_start = 100
        # read covers genomic [104, 107) → region indices 4,5,6
        rstart = 104
        rseq = "TGA"  # length 3

        scores, covered = score_read_against_region(
            rstart=rstart,
            rseq=rseq,
            region_start=region_start,
            c_type="allC",
            no_endog=False,
            ref=ref,
        )
        self.assertIsNotNone(scores)
        self.assertEqual(len(scores), len(ref))
        self.assertEqual(len(covered), len(ref))

        expected_covered = np.array(
            [False, False, False, False, True, True, True, False],
            dtype=bool,
        )
        np.testing.assert_array_equal(covered, expected_covered)

    def test_CG_mode_NCG_scoring(self):
        """
        c_type='CG' → score the C in 'CG' (NCG motif conceptually).
        """
        ref = "ACGTCGTA"
        # indices: 0 A,1 C,2 G,3 T,4 C,5 G,6 T,7 A
        # CG at (1,2) and (4,5) → score positions 1 and 4
        scores, covered = score_read_against_region(
            rstart=0,
            rseq=ref,          # perfectly matches reference
            region_start=0,
            c_type="CG",
            no_endog=False,
            ref=ref,
        )
        self.assertIsNotNone(scores)
        expected = np.array(
            [-1, 0, -1, -1, 0, -1, -1, -1],
            dtype=np.int8,
        )
        np.testing.assert_array_equal(scores, expected)
        self.assertTrue(covered.all())

        # If we mutate the C at position 4 to T, that position should become 1
        rseq_mut = list(ref)
        rseq_mut[4] = "T"
        rseq_mut = "".join(rseq_mut)
        scores_mut, _ = score_read_against_region(
            rstart=0,
            rseq=rseq_mut,
            region_start=0,
            c_type="CG",
            no_endog=False,
            ref=ref,
        )
        expected_mut = expected.copy()
        expected_mut[4] = 1
        np.testing.assert_array_equal(scores_mut, expected_mut)

    def test_GC_mode_GCN_no_endog_true(self):
        """
        c_type='GC', no_endog=True → GCN (including GCG).
        Score the C of the GC* triplet.
        """
        ref = "AGCGCAGCG"
        # indices: 0 A,1 G,2 C,3 G,4 C,5 A,6 G,7 C,8 G
        # GC at (1,2) triplet GCG → C at 2
        # GC at (3,4) triplet GCA → C at 4
        # GC at (6,7) triplet GCG → C at 7
        scores, covered = score_read_against_region(
            rstart=0,
            rseq=ref,
            region_start=0,
            c_type="GC",
            no_endog=True,     # include both GCG and GCA contexts
            ref=ref,
        )
        expected = np.array(
            [-1, -1, 0, -1, 0, -1, -1, 0, -1],
            dtype=np.int8,
        )
        np.testing.assert_array_equal(scores, expected)
        self.assertTrue(covered.all())

    def test_GC_mode_GCH_no_endog_false(self):
        """
        c_type='GC', no_endog=False → GCH (exclude GCG).
        So only C where the following base is not G.
        """
        ref = "AGCGCAGCG"
        # same as above; only the GCA at indices 3–5 (C at 4) should be scored
        # The GCGs at 1–3 and 6–8 should be excluded.
        # Make read all 'C' at those positions to test scoring.
        rseq = ref
        scores, _ = score_read_against_region(
            rstart=0,
            rseq=rseq,
            region_start=0,
            c_type="GC",
            no_endog=False,   # exclude GCG
            ref=ref,
        )
        expected = np.array(
            [-1, -1, -1, -1, 0, -1, -1, -1, -1],
            dtype=np.int8,
        )
        np.testing.assert_array_equal(scores, expected)

        # If we change the C at position 4 to A, that site should become 1
        rseq_mut = list(ref)
        rseq_mut[4] = "A"
        rseq_mut = "".join(rseq_mut)
        scores_mut, _ = score_read_against_region(
            rstart=0,
            rseq=rseq_mut,
            region_start=0,
            c_type="GC",
            no_endog=False,
            ref=ref,
        )
        expected_mut = expected.copy()
        expected_mut[4] = 1
        np.testing.assert_array_equal(scores_mut, expected_mut)

    def test_both_dimers_union_of_CG_and_GC(self):
        """
        c_type='both_dimers' → union of CG- and GC-derived positions.
        """
        ref = "AGCCG"
        # indices: 0 A, 1 G, 2 C, 3 C, 4 G
        # Union of C positions = {2,3}
        scores, _ = score_read_against_region(
            rstart=0,
            rseq=ref,
            region_start=0,
            c_type="both_dimers",
            no_endog=True,  # GC side → GCN (doesn't matter here)
            ref=ref,
        )
        expected = np.array(
            [-1, -1, 0, 0, -1],
            dtype=np.int8,
        )
        np.testing.assert_array_equal(scores, expected)

    def test_allC_scores_every_cytosine(self):
        """c_type='allC' → every C in the region is scored."""
        ref = "ACCTAGC"
        # indices: 0 A,1 C,2 C,3 T,4 A,5 G,6 C
        rseq = ref
        scores, covered = score_read_against_region(
            rstart=0,
            rseq=rseq,
            region_start=0,
            c_type="allC",
            no_endog=False,
            ref=ref,
        )
        expected = np.array(
            [-1, 0, 0, -1, -1, -1, 0],
            dtype=np.int8,
        )
        np.testing.assert_array_equal(scores, expected)
        self.assertTrue(covered.all())

        # Change one C → non-C and check it flips to 1
        rseq_mut = list(ref)
        rseq_mut[2] = "T"
        rseq_mut = "".join(rseq_mut)
        scores_mut, _ = score_read_against_region(
            rstart=0,
            rseq=rseq_mut,
            region_start=0,
            c_type="allC",
            no_endog=False,
            ref=ref,
        )
        expected_mut = expected.copy()
        expected_mut[2] = 1
        np.testing.assert_array_equal(scores_mut, expected_mut)


if __name__ == "__main__":
    unittest.main()