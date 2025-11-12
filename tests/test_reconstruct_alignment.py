# test_reconstruct_alignment.py
import os
import tempfile
import unittest
from typing import Tuple, Dict
import pysam
from dSMF_footprints_clustering_py3 import reconstruct_aligned_query


def _write_bam(path: str) -> None:
    """
    Create a tiny BAM (and .bai) with several CIGAR cases:
      rM:   10M
      rSDM: 2S 5M 3D 4M 2S
      rMI:  5M 2I 5M
      rMix: 3M 1I 2M 2D 4M
    All on chr1, which we declare length 1000.
    """
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        # rM: 10M @ pos 100
        a = pysam.AlignedSegment()
        a.query_name = "rM"
        a.query_sequence = "ACGTACGTAC"  # 10 bases
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 100
        a.cigar = [(0, 10)]              # 10M
        a.mapping_quality = 60
        bam.write(a)

        # rSDM: 2S 5M 3D 4M 2S @ pos 200
        # soft clips are consumed from query but not written; deletions -> 'N'*len
        b = pysam.AlignedSegment()
        b.query_name = "rSDM"
        b.query_sequence = "TT" + "GGGGG" + "CCCC" + "AA"  # 2S + 5M + 4M + 2S
        b.flag = 0
        b.reference_id = 0
        b.reference_start = 200
        b.cigar = [(4, 2), (0, 5), (2, 3), (0, 4), (4, 2)]  # S M D M S
        b.mapping_quality = 60
        bam.write(b)

        # rMI: 5M 2I 5M @ pos 300
        c = pysam.AlignedSegment()
        c.query_name = "rMI"
        c.query_sequence = "AAAAA" + "GG" + "TTTTT"  # M + I + M
        c.flag = 0
        c.reference_id = 0
        c.reference_start = 300
        c.cigar = [(0, 5), (1, 2), (0, 5)]
        c.mapping_quality = 60
        bam.write(c)

        # rMix: 3M 1I 2M 2D 4M @ pos 400
        d = pysam.AlignedSegment()
        d.query_name = "rMix"
        d.query_sequence = "TTT" + "G" + "AA" + "CCCC"
        d.flag = 0
        d.reference_id = 0
        d.reference_start = 400
        d.cigar = [(0, 3), (1, 1), (0, 2), (2, 2), (0, 4)]
        d.mapping_quality = 60
        bam.write(d)

    pysam.index(path)


class TestReconstructAlignedQueryIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.bam_path = os.path.join(cls.tmpdir.name, "mini.bam")
        _write_bam(cls.bam_path)
        cls.bam = pysam.AlignmentFile(cls.bam_path, "rb")

        # Expected reconstructed strings:
        # - S (soft clip): skipped
        # - I (insertion): skipped
        # - D (deletion): pad with 'N'
        cls.expected: Dict[str, Tuple[int, str]] = {
            "rM":   (100, "ACGTACGTAC"),             # 10M -> copy sequence
            "rSDM": (200, "GGGGG" + "NNN" + "CCCC"), # 5M + 3D + 4M
            "rMI":  (300, "AAAAA" + "TTTTT"),        # insertion removed
            "rMix": (400, "TTT" + "AA" + "NN" + "CCCC"),
        }

    @classmethod
    def tearDownClass(cls):
        cls.bam.close()
        # clean up BAM and index
        bai = cls.bam_path + ".bai"
        if os.path.exists(bai):
            os.remove(bai)
        cls.tmpdir.cleanup()

    def test_each_alignment_matches_expected(self):
        seen = set()
        for aln in self.bam.fetch(until_eof=True):
            qn = aln.query_name
            self.assertIn(qn, self.expected, f"Unexpected read {qn} in test BAM")
            rstart, seq = reconstruct_aligned_query(aln)

            exp_start, exp_seq = self.expected[qn]
            self.assertEqual(rstart, exp_start, f"{qn}: start mismatch")
            self.assertEqual(seq, exp_seq, f"{qn}: reconstructed seq mismatch")

            # Also validate length == ref span
            ref_span = aln.reference_end - aln.reference_start
            self.assertEqual(len(seq), ref_span, f"{qn}: length != ref span")

            seen.add(qn)

        # Ensure we tested every synthetic read we wrote
        self.assertEqual(seen, set(self.expected.keys()))


if __name__ == "__main__":
    unittest.main()
