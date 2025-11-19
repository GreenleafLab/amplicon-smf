# tests/test_merge_alignments.py
import unittest
import numpy as np
from dSMF_footprints_clustering_py3 import merge_alignments

class TestMergeAlignments(unittest.TestCase):

    def setUp(self):
        # Two forward segments from the SAME template with a 4bp overlap.
        # rstart: 0 (first segment), 196 (second segment)
        self.seq1 = (
            "GGTGAGTAGTAGGGTAGAAGATGGTATTGTGGAAGTTATGTATGTTTGCTGCAATTGGATTTAGCTATGTTAGGCTTTGCAGTGTCTTAGTTATTTGAGTT"
            "TGGTGTTATGGTATTGTAGTTGTTAGTGGGTTTATTTTGTTAATAGAGTTTTGTTAGTTTGTAAGTATGTTATATAGGTTATTTGGGTAGCTATATGTA"
        )
        self.seq2 = (
            "TGTATGTTGAGGTATGATATGTTTTTTTTAGTAAGCTATTGTTGTTAGGATGCAGGGAGGCTTTTGCTTTTGAGTGCATTGCTTATGGTAATTTTAGTAT"
            "AGTGTGTTTGGAGAGTTTTGTATTAGTTGGAGTATTGGTTGTATAAATTGTTGTTAAGGTAATTGGAATGGGTNAATTGGGTGGAAGGATGTAGAGTGATG"
        )
        # sanity: lengths are ~200
        self.len1 = len(self.seq1)
        self.len2 = len(self.seq2)
        self.rstart1 = 0
        self.rstart2 = 196  # 4 bp overlap if each is length 200

        # Region starts at 0; reference needs to be at least as long as the merged span
        max_end = max(self.rstart1 + self.len1, self.rstart2 + self.len2)
        self.region_start = 0
        # Reference here is arbitrary (no conflicts test), but long enough
        self.ref = "A" * max_end

    def test_simple_two_segment_merge_no_conflict(self):
        segs = [(self.rstart1, self.seq1), (self.rstart2, self.seq2)]
        rstart_m, seq_m, cov = merge_alignments(segs, self.region_start, self.ref, prefer_ref_on_conflict=True)

        # merged starts where the first begins
        self.assertEqual(rstart_m, 0)

        # Length should cover from 0 to max_end (no trimming at ends since both start/end in non-Ns)
        expected_len = max(self.rstart1 + self.len1, self.rstart2 + self.len2) - rstart_m
        self.assertEqual(len(seq_m), expected_len)
        self.assertEqual(len(cov), expected_len)

        # Coverage should be >=1 across the union; positions in the 4bp overlap are 2
        # Compute “overlap positions with both non-N”
        o_len = max(0, (self.rstart1 + self.len1) - self.rstart2)
        self.assertGreaterEqual(o_len, 0)
        # Slices of the overlapping part (relative to segment coordinates)
        ovl1 = self.seq1[self.len1 - o_len : self.len1] if o_len > 0 else ""
        ovl2 = self.seq2[:o_len] if o_len > 0 else ""
        nonN_both = sum(1 for a, b in zip(ovl1, ovl2) if a != "N" and b != "N")
        # positions of overlap in merged coordinates: [rstart2, rstart2+o_len)
        two_cov = int(np.sum(cov[self.rstart2 : self.rstart2 + o_len] == 2))
        self.assertEqual(two_cov, nonN_both)

        # Sanity: merged contains the left chunk exactly at the left side
        self.assertTrue(seq_m.startswith(self.seq1[:50]))  # sample prefix check

    def test_conflict_resolution_pref_ref(self):
        # Create a tiny synthetic conflict at the overlap:
        # seq1 tail ends with 'A' at overlap pos; seq2 head is 'C' there.
        s1 = "TTTTTAAAA"        # length 9
        s2 = "CTTTTTTTT"        # length 9; conflict at merged index rstart2 = 5, we overlap 4 bp
        r1, r2 = 0, 5
        ref = "TTTTTCTTTTT"     # at position 5, ref base = 'C' (matches s2)
        segs = [(r1, s1), (r2, s2)]
        rstart_m, seq_m, cov = merge_alignments(segs, 0, ref, prefer_ref_on_conflict=True)

        # At conflict index (5), prefer ref -> picks 'C' (s2)
        self.assertEqual(seq_m[5], "C")
        # Coverage at the first overlap position should be 2 (both informative)
        self.assertEqual(cov[5], 2)

    def test_conflict_resolution_no_prefer(self):
        s1 = "TTTTTAAAA"
        s2 = "CTTTTTTTT"
        r1, r2 = 0, 5
        ref = "TTTTTATTTTT"     # ref pos5 = 'A' (matches s1), but we won’t prefer
        segs = [(r1, s1), (r2, s2)]
        rstart_m, seq_m, cov = merge_alignments(segs, 0, ref, prefer_ref_on_conflict=False)

        # With no preference, conflicting informative bases become 'N'
        self.assertEqual(seq_m[5], "N")
        self.assertEqual(cov[5], 2)  # coverage still 2

    def test_all_N_trims_to_empty(self):
        segs = [(0, "NNNNN"), (3, "NNN")]
        rstart_m, seq_m, cov = merge_alignments(segs, 0, "A" * 10, prefer_ref_on_conflict=True)
        self.assertEqual(seq_m, "")
        self.assertEqual(len(cov), 0)
        self.assertEqual(rstart_m, 0)  # returns region_start

    def test_leading_trailing_N_trim(self):
        segs = [(0, "NNNACTGNNN")]  # inside there is ACTG surrounded by Ns
        rstart_m, seq_m, cov = merge_alignments(segs, 0, "A" * 20)
        self.assertEqual(seq_m, "ACTG")
        self.assertEqual(rstart_m, 3)
        self.assertTrue(np.all(cov == 1))

if __name__ == "__main__":
    unittest.main()